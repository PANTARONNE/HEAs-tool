#!/bin/bash
#SBATCH -N 1
#SBATCH --ntasks-per-node=4
#SBATCH -c 8 #此处请不要改动
#SBATCH -J
#SBATCH --gres=dcu:4
#SBATCH -p wzhdnormal
#SBATCH --output=out_%j.log
#SBATCH --error=err_%j.log

set -uo pipefail

# --- self-resubmitting relaxation chain --------------------------------------
# After VASP finishes, check-convergence.sh classifies the run. If it merely ran
# out of ionic steps we archive this generation, roll CONTCAR into POSCAR and
# resubmit ourselves, up to MAX_GEN generations. A crash or an ill-behaved SCF
# stops the chain and leaves a marker file for the driver to notice.
MAX_GEN=${MAX_GEN:-2}
workdir=${SLURM_SUBMIT_DIR:-$PWD}
cd "$workdir" || exit 1
check_script="${workdir}/check-convergence.sh"

for marker in CONVERGED NOT_CONVERGED_MAX_GEN CRASHED SCF_ILL CHECK_ERROR; do
    if [[ -e "$marker" ]]; then
        echo "[chain] terminal marker ${marker} already exists; exiting without running VASP." >&2
        exit 0
    fi
done

module purge
source /work/home/acsppswc7s/apprepo/vasp/642-dtk23.10_all/scripts/env.sh

#指定版本,std,gam或ncl版
VASP_EXE=vasp_gam

config=config.${SLURM_JOB_ID}
echo -e "-genv OMP_NUM_THREADS 6 \c" > $config
for i in $(scontrol show hostnames $SLURM_NODELIST)
do
  for ((j=0; j<4; j++))
  do
    echo "-host $i -env HIP_VISIBLE_DEVICES $j -n 1 numactl --cpunodebind=$j --membind=$j $VASP_EXE" >> $config
  done
done

ulimit -s unlimited

export NCCL_IB_HCA="mlx5_0"
export HSA_FORCE_FINE_GRAIN_PCIE=1

mpirun -configfile $config

# --- post-run classification and chaining ------------------------------------
gen=$(cat .gen_count 2>/dev/null || echo 1)

# Archive this generation's key outputs so the resubmit chain never overwrites
# earlier history.
gen_dir=$(printf 'gen_%02d' "$gen")
mkdir -p "$gen_dir"
for f in OUTCAR OSZICAR CONTCAR XDATCAR vasprun.xml "out_${SLURM_JOB_ID}.log" "err_${SLURM_JOB_ID}.log"; do
    [[ -e "$f" ]] && cp -p -- "$f" "$gen_dir/" 2>/dev/null || true
done

if [[ ! -x "$check_script" && ! -f "$check_script" ]]; then
    echo "[chain] check-convergence.sh not found at $check_script" >&2
    echo "check-convergence.sh missing" > CHECK_ERROR
    exit 1
fi

bash "$check_script" "$workdir"
status=$?

case "$status" in
    0)
        echo "[chain] gen ${gen}: converged." >&2
        touch CONVERGED
        ;;
    10)
        if [[ "$gen" -ge "$MAX_GEN" ]]; then
            echo "[chain] gen ${gen}: still not converged at MAX_GEN=${MAX_GEN}; stopping." >&2
            touch NOT_CONVERGED_MAX_GEN
        elif [[ ! -s CONTCAR ]]; then
            echo "[chain] gen ${gen}: CONTCAR missing/empty, cannot continue." >&2
            touch CRASHED
        else
            resubmit_marker=".resubmitted_gen_${gen}"
            if [[ -e "$resubmit_marker" ]]; then
                echo "[chain] gen ${gen}: resubmission was already requested; not submitting again." >&2
                exit 0
            fi
            next=$((gen + 1))
            echo "[chain] gen ${gen}: not converged; resubmitting as gen ${next}." >&2
            cp -p CONTCAR POSCAR
            echo "$next" > .gen_count
            touch "$resubmit_marker"
            if ! sbatch --export=ALL,MAX_GEN="$MAX_GEN" "$0"; then
                rm -f "$resubmit_marker"
                echo "[chain] gen ${gen}: sbatch failed." >&2
                touch CHECK_ERROR
                exit 1
            fi
        fi
        ;;
    20)
        echo "[chain] gen ${gen}: VASP crashed/timed out; not resubmitting." >&2
        touch CRASHED
        ;;
    30)
        echo "[chain] gen ${gen}: SCF ill-behaved; not resubmitting." >&2
        touch SCF_ILL
        ;;
    *)
        echo "[chain] gen ${gen}: unexpected check status ${status}." >&2
        touch CHECK_ERROR
        ;;
esac

rm -f "$config"
