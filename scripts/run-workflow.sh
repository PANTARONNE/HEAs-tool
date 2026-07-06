#!/usr/bin/env bash
#
# run-workflow.sh — end-to-end driver for the HEA DFT dataset workflow described
# in CALCULATION_WORKFLOW.md. It runs on the Slurm login node and blocks until
# each stage finishes, driving:
#
#   1. create-sample        -> register a new HEA composition (initial slab)
#   2. slab-relax (VASP)     -> self-resubmitting relaxation, CONTCAR -> cif,
#                               record-relaxed into the dataset
#   3. index-surface         -> grid the surface atoms
#   4. hamilton (OpenMX)     -> SCF, then extract-hamiltonian into the dataset
#
# Usage:
#   run-workflow.sh -e Fe Co Ni Cr Mn [options]
#   run-workflow.sh --surface-id Co_13-Cr_13-Fe_13-Mn_12-Ni_13 [options]   # skip step 1
#
# Options:
#   -e, --elements EL...     elements for create-sample (step 1)
#   -r, --ratios R...        optional ratios (random when omitted)
#       --surface-id ID      use an existing composition; skip create-sample
#       --root DIR           dataset root (default: dataset)
#       --workspace DIR      calculation workspace (default: workspace)
#       --openmx-data DIR    server-side OpenMX DFT_DATA path (default: DFT_DATA19)
#       --poll SECONDS       job-state poll interval (default: 60)
#       --max-gen N          max VASP relaxation generations (default: 5)
#       --skip-hamilton      stop after step 3
#   -h, --help

set -Eeuo pipefail

script_dir=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)
repo_dir=$(cd -- "${script_dir}/.." && pwd)
PYTHON=${PYTHON:-python}

# --- defaults ----------------------------------------------------------------
root="dataset"
workspace="workspace"
openmx_data="DFT_DATA19"
poll=60
max_gen=5
skip_hamilton=0
surface_id=""
elements=()
ratios=()

die() { printf 'Error: %s\n' "$*" >&2; exit 1; }
log() { printf '[workflow] %s\n' "$*" >&2; }

# --- argument parsing --------------------------------------------------------
while [[ $# -gt 0 ]]; do
    case "$1" in
        -e|--elements)
            shift
            while [[ $# -gt 0 && "$1" != -* ]]; do elements+=("$1"); shift; done
            ;;
        -r|--ratios)
            shift
            while [[ $# -gt 0 && "$1" != -* ]]; do ratios+=("$1"); shift; done
            ;;
        --surface-id) surface_id=$2; shift 2 ;;
        --root)       root=$2; shift 2 ;;
        --workspace)  workspace=$2; shift 2 ;;
        --openmx-data) openmx_data=$2; shift 2 ;;
        --poll)       poll=$2; shift 2 ;;
        --max-gen)    max_gen=$2; shift 2 ;;
        --skip-hamilton) skip_hamilton=1; shift ;;
        -h|--help)    grep '^#' "$0" | sed 's/^# \{0,1\}//'; exit 0 ;;
        *) die "unknown argument: $1" ;;
    esac
done

command -v sbatch >/dev/null 2>&1 || die "sbatch not found; run on a Slurm login node."
[[ -f "${repo_dir}/hea_dataset.py" ]] || die "hea_dataset.py not found in ${repo_dir}."

# Wait until a Slurm job leaves the queue (completed, failed or cancelled).
wait_for_job() {
    local jid=$1
    log "waiting for Slurm job ${jid} ..."
    while squeue -h -j "$jid" 2>/dev/null | grep -q .; do
        sleep "$poll"
    done
}

# Wait until one of the given marker files appears in a directory.
wait_for_marker() {
    local dir=$1; shift
    log "waiting for completion markers in ${dir} ..."
    while true; do
        for m in "$@"; do
            if [[ -e "${dir}/${m}" ]]; then
                printf '%s\n' "$m"
                return 0
            fi
        done
        sleep "$poll"
    done
}

# ============================================================================
# Step 1: create and register the initial HEA structure
# ============================================================================
if [[ -z "$surface_id" ]]; then
    [[ ${#elements[@]} -gt 0 ]] || die "provide --elements or --surface-id."
    log "step 1: create-sample -e ${elements[*]} ${ratios:+-r ${ratios[*]}}"
    create_args=(create-sample --root "$root" -e "${elements[@]}")
    [[ ${#ratios[@]} -gt 0 ]] && create_args+=(-r "${ratios[@]}")
    create_out=$("$PYTHON" "${repo_dir}/hea_dataset.py" "${create_args[@]}")
    printf '%s\n' "$create_out" >&2
    surface_id=$(printf '%s\n' "$create_out" \
        | awk -F':' '/Surface ID/ { gsub(/[[:space:]]/, "", $2); print $2; exit }')
    [[ -n "$surface_id" ]] || die "could not parse Surface ID from create-sample output."
else
    log "step 1: skipped (using existing surface-id ${surface_id})."
fi
log "surface-id = ${surface_id}"

sample_dir="${workspace}/${surface_id}"
initial_cif="${root}/${surface_id}/structures/00_initial_sqs.cif"
[[ -f "$initial_cif" ]] || die "initial structure not found: ${initial_cif}"

# ============================================================================
# Step 2: relax the slab with VASP, then record the relaxed structure
# ============================================================================
relax_dir="${sample_dir}/slab-relax"
log "step 2: VASP slab relaxation in ${relax_dir}"
mkdir -p "$relax_dir"
cp -- "$initial_cif" "${relax_dir}/${surface_id}.cif"
cp -- "${script_dir}/check-convergence.sh" "${relax_dir}/check-convergence.sh"
chmod +x "${relax_dir}/check-convergence.sh"

(
    cd "$relax_dir"
    # vasp-inputs.sh builds INCAR/POSCAR/POTCAR/KPOINTS and a job-named slurm
    # script; the surface-id has no spaces so it is a valid Slurm job name.
    bash "${script_dir}/vasp-inputs.sh" "${surface_id}.cif"
    echo 1 > .gen_count
    rm -f CONVERGED NOT_CONVERGED_MAX_GEN CRASHED SCF_ILL CHECK_ERROR
    MAX_GEN="$max_gen" sbatch vasp-gam.slurm
)

marker=$(wait_for_marker "$relax_dir" \
    CONVERGED NOT_CONVERGED_MAX_GEN CRASHED SCF_ILL CHECK_ERROR)
case "$marker" in
    CONVERGED) log "step 2: relaxation converged." ;;
    *) die "step 2: relaxation did not converge (marker: ${marker}). See ${relax_dir}." ;;
esac

# Convert the final CONTCAR back to CIF. VASPKIT/VASP keep POSCAR's per-element
# grouping through CONTCAR, so ASE reproduces the atom order record-relaxed
# expects (same count and per-atom element sequence as the initial slab).
relaxed_cif="${relax_dir}/01_relaxed_slab.cif"
log "step 2: converting CONTCAR -> ${relaxed_cif}"
"$PYTHON" - "$relax_dir/CONTCAR" "$relaxed_cif" <<'PY'
import sys
from ase.io import read, write
write(sys.argv[2], read(sys.argv[1], format="vasp"))
PY

log "step 2: record-relaxed"
"$PYTHON" "${repo_dir}/hea_dataset.py" record-relaxed \
    --root "$root" --surface-id "$surface_id" --relaxed-cif "$relaxed_cif" >&2

# ============================================================================
# Step 3: index the surface atoms
# ============================================================================
log "step 3: index-surface"
"$PYTHON" "${repo_dir}/hea_dataset.py" index-surface \
    --root "$root" --surface-id "$surface_id" >&2

if [[ "$skip_hamilton" -eq 1 ]]; then
    log "done (stopped after step 3 as requested). surface-id=${surface_id}"
    exit 0
fi

# ============================================================================
# Step 4: OpenMX Hamiltonian calculation, then extract-hamiltonian
# ============================================================================
hamilton_dir="${sample_dir}/hamilton"
job_name="${surface_id}-hamilton"
relaxed_registered="${root}/${surface_id}/structures/01_relaxed_slab.cif"
log "step 4: OpenMX Hamiltonian in ${hamilton_dir}"
mkdir -p "$hamilton_dir"

# input.dat -> System.Name "input" -> input.scfout (matches extract-hamiltonian).
"$PYTHON" "${repo_dir}/cif_to_openmx.py" "$relaxed_registered" \
    --data-path "$openmx_data" -o "${hamilton_dir}/input.dat" >&2

# Copy the OpenMX slurm template, set the job name and the input .dat filename.
awk -v job="$job_name" '
    { sub(/\r$/, "") }
    /^[[:space:]]*#SBATCH[[:space:]]+-J/ { print "#SBATCH -J " job; next }
    /openmx[[:space:]]+[^[:space:]]+\.dat/ { sub(/openmx[[:space:]]+[^[:space:]]+\.dat/, "openmx input.dat") }
    { print }
' "${script_dir}/openmx.slurm" > "${hamilton_dir}/openmx.slurm"

hjid=$(
    cd "$hamilton_dir"
    sbatch --parsable openmx.slurm
)
[[ -n "$hjid" ]] || die "step 4: failed to submit OpenMX job."
wait_for_job "$hjid"

scfout="${hamilton_dir}/input.scfout"
[[ -s "$scfout" ]] || die "step 4: OpenMX did not produce ${scfout}. See ${hamilton_dir}."

log "step 4: extract-hamiltonian"
"$PYTHON" "${repo_dir}/hea_dataset.py" extract-hamiltonian \
    --root "$root" --surface-id "$surface_id" \
    --dat "${hamilton_dir}/input.dat" --scfout "$scfout" >&2

log "done. surface-id=${surface_id}"
