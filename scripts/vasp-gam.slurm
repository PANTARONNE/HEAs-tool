#!/bin/bash
#SBATCH -N 1
#SBATCH --ntasks-per-node=4
#SBATCH -c 8 #此处请不要改动
#SBATCH -J
#SBATCH --gres=dcu:4
#SBATCH -p wzhdnormal
#SBATCH --output=out_%j.log
#SBATCH --error=err_%j.log

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
