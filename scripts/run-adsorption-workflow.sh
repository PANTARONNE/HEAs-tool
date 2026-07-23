#!/usr/bin/env bash
#
# run-adsorption-workflow.sh — end-to-end driver for the adsorption-energy
# workflow described in docs/CALCULATION_WORKFLOW_ADSORPTION.md. It runs on the Slurm login
# node, on top of a surface that already completed the Hamiltonian workflow
# (slab relaxed + total energy recorded), and drives:
#
#   1. detect-sites             -> enumerate FCC hollow sites on the relaxed slab
#   2. create-adsorbate-records -> build one initial adsorbate structure per site
#   3. slab+ads relax (VASP)    -> self-resubmitting relaxation per site,
#                                  CONTCAR -> cif
#   4. adsorption energy         -> E_ads = E(slab+ads) - E(slab) - E(ref),
#                                  then record-energy into the dataset
#
# Usage:
#   run-adsorption-workflow.sh --surface-id ID -a N [options]
#   run-adsorption-workflow.sh --surface-id ID -a NH2 --site 5 [options]
#
# Options:
#   --surface-id ID       existing composition that finished workflow 1 (required)
#   -a, --adsorbate SP    adsorbate species: N | NH | NH2 | NH3 (required)
#   --site SEL...         site numbers to run (e.g. --site 1 3 5), or omit / "all"
#                         for every detected site (default: all)
#       --root DIR        dataset root (default: dataset)
#       --workspace DIR   calculation workspace (default: workspace)
#       --slab-energy E   override E(slab) in eV (default: read from the dataset)
#       --ref-energy E    override the adsorbate reference energy in eV for this run
#       --height A        adsorbate height above the hollow plane (default: 1.25)
#       --nh A            N-H bond length for NHx (default: 1.02)
#       --poll SECONDS    job-state poll interval (default: 60)
#       --max-gen N       max VASP relaxation generations per site (default: 5)
#       --skip-detect     reuse an existing fcc_sites.jsonl (skip step 1)
#       --force           resubmit sites that already carry a terminal marker
#   -h, --help
#
# ============================================================================
# ADSORBATE REFERENCE ENERGIES  --  FILL THESE IN BEFORE USE
# ----------------------------------------------------------------------------
# Adsorption energy is defined as
#     E_ads = E(slab + adsorbate)  -  E(slab)  -  E_ref(adsorbate)
# where E_ref is the reference (e.g. gas-phase) total energy of the adsorbate
# fragment in eV. Put the correct value for each species you intend to use;
# leave the others empty. A run aborts if the species it needs is still empty
# (unless --ref-energy is supplied on the command line).
# ============================================================================
declare -A ADSORBATE_REF_ENERGY=(
    [N]="-8.3177"      # eV  e.g. 0.5 * E(N2)
    [NH]="-8.1023"     # eV
    [NH2]="-13.5359"    # eV
    [NH3]="-19.5405"    # eV
)

set -Eeuo pipefail

script_dir=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)
repo_dir=$(cd -- "${script_dir}/.." && pwd)
PYTHON=${PYTHON:-python}

# --- defaults ----------------------------------------------------------------
root="dataset"
workspace="workspace"
surface_id=""
adsorbate=""
sites=()
slab_energy=""
ref_energy=""
height=1.25
nh=1.02
poll=60
max_gen=5
skip_detect=0
force=0

die() { printf 'Error: %s\n' "$*" >&2; exit 1; }
log() { printf '[ads-workflow] %s\n' "$*" >&2; }

# True while a Slurm job with the given name is still queued/running. The VASP
# self-resubmit chain keeps the same job name across generations, so this stays
# true for the whole chain and only turns false once the chain has left the
# queue (finished, failed, cancelled or killed).
job_alive() {
    squeue -h -u "${USER:-$(id -un)}" -n "$1" 2>/dev/null | grep -q .
}

# --- argument parsing --------------------------------------------------------
while [[ $# -gt 0 ]]; do
    case "$1" in
        --surface-id) surface_id=$2; shift 2 ;;
        -a|--adsorbate) adsorbate=$2; shift 2 ;;
        --site)
            shift
            while [[ $# -gt 0 && "$1" != -* ]]; do sites+=("$1"); shift; done
            ;;
        --root)        root=$2; shift 2 ;;
        --workspace)   workspace=$2; shift 2 ;;
        --slab-energy) slab_energy=$2; shift 2 ;;
        --ref-energy)  ref_energy=$2; shift 2 ;;
        --height)      height=$2; shift 2 ;;
        --nh)          nh=$2; shift 2 ;;
        --poll)        poll=$2; shift 2 ;;
        --max-gen)     max_gen=$2; shift 2 ;;
        --skip-detect) skip_detect=1; shift ;;
        --force)       force=1; shift ;;
        -h|--help)     awk 'NR==1{next} /^#/{sub(/^# ?/,""); print; next} {exit}' "$0"; exit 0 ;;
        *) die "unknown argument: $1" ;;
    esac
done

[[ -n "$surface_id" ]] || die "provide --surface-id."
[[ -n "$adsorbate" ]] || die "provide -a/--adsorbate (N|NH|NH2|NH3)."
case "$adsorbate" in
    N|NH|NH2|NH3) ;;
    *) die "unsupported adsorbate: ${adsorbate} (choose N, NH, NH2 or NH3)." ;;
esac

command -v sbatch >/dev/null 2>&1 || die "sbatch not found; run on a Slurm login node."
[[ -f "${repo_dir}/hea_dataset.py" ]] || die "hea_dataset.py not found in ${repo_dir}."

relaxed_registered="${root}/${surface_id}/structures/01_relaxed_slab.cif"
[[ -f "$relaxed_registered" ]] || die \
    "relaxed slab not found: ${relaxed_registered}. Run workflow 1 (run-hamiltonian-workflow.sh) first."

# --- resolve the adsorbate reference energy ----------------------------------
if [[ -z "$ref_energy" ]]; then
    ref_energy=${ADSORBATE_REF_ENERGY[$adsorbate]:-}
fi
[[ -n "$ref_energy" ]] || die \
    "no reference energy for ${adsorbate}. Fill ADSORBATE_REF_ENERGY[${adsorbate}] near the top of this script, or pass --ref-energy."

# --- resolve E(slab) ---------------------------------------------------------
if [[ -z "$slab_energy" ]]; then
    log "reading E(slab) from the dataset index ..."
    slab_energy=$(
        "$PYTHON" - "$root" "$surface_id" <<'PY'
import sqlite3, sys
from pathlib import Path
root, surface_id = sys.argv[1], sys.argv[2]
db = Path(root) / "index.sqlite"
value = ""
if db.is_file():
    with sqlite3.connect(db) as con:
        row = con.execute(
            "SELECT total_energy_eV FROM surfaces WHERE surface_id=?",
            (surface_id,),
        ).fetchone()
        if row and row[0] is not None:
            value = repr(float(row[0]))
print(value)
PY
    )
fi
[[ -n "$slab_energy" ]] || die \
    "E(slab) is not recorded for ${surface_id}. Re-run workflow 1 record-relaxed with --energy, or pass --slab-energy."
log "E(slab)      = ${slab_energy} eV"
log "E(ref) ${adsorbate}   = ${ref_energy} eV"

# ============================================================================
# Step 1: detect FCC hollow sites on the relaxed slab
# ============================================================================
fcc_sites="${root}/${surface_id}/metadata/fcc_sites.jsonl"
if [[ "$skip_detect" -eq 1 && -f "$fcc_sites" ]]; then
    log "step 1: skipped (reusing ${fcc_sites})."
else
    log "step 1: detect-sites"
    "$PYTHON" "${repo_dir}/hea_dataset.py" detect-sites \
        --root "$root" --surface-id "$surface_id" >&2
fi
[[ -f "$fcc_sites" ]] || die "no FCC sites file produced: ${fcc_sites}"

# ============================================================================
# Step 2: build one initial adsorbate structure per site
# ============================================================================
log "step 2: create-adsorbate-records -a ${adsorbate}"
"$PYTHON" "${repo_dir}/hea_dataset.py" create-adsorbate-records \
    --root "$root" --surface-id "$surface_id" \
    --adsorbate "$adsorbate" --height "$height" --nh "$nh" >&2

# --- resolve which sites to process ------------------------------------------
all_sites=()
while IFS= read -r n; do
    [[ -n "$n" ]] && all_sites+=("$n")
done < <(
    "$PYTHON" - "$fcc_sites" <<'PY'
import json, sys
idx = []
with open(sys.argv[1]) as fh:
    for line in fh:
        line = line.strip()
        if line:
            idx.append(int(json.loads(line)["site_index"]))
for i in sorted(idx):
    print(i)
PY
)
[[ ${#all_sites[@]} -gt 0 ]] || die "no site indices parsed from ${fcc_sites}"

if [[ ${#sites[@]} -eq 0 || "${sites[0]}" == "all" ]]; then
    sites=("${all_sites[@]}")
else
    for s in "${sites[@]}"; do
        [[ "$s" =~ ^[0-9]+$ ]] || die "invalid --site value: ${s}"
        printf '%s\n' "${all_sites[@]}" | grep -qx "$s" \
            || die "site ${s} is not in the detected range (1..${#all_sites[@]})."
    done
fi
log "processing ${#sites[@]} site(s): ${sites[*]}"

# ============================================================================
# Step 3: relax each slab+adsorbate structure with VASP (parallel submit)
# ============================================================================
declare -a site_dirs
declare -a site_jobs
ads_ws="${workspace}/${surface_id}/adsorbate/${adsorbate}"

site_id_of() { printf 'site_%04d' "$1"; }

submitted=0
for s in "${sites[@]}"; do
    site_id=$(site_id_of "$s")
    init_cif="${root}/${surface_id}/adsorbates/${adsorbate}/${site_id}/00_initial_adsorbate.cif"
    [[ -f "$init_cif" ]] || die "initial adsorbate structure missing: ${init_cif}"

    dir="${ads_ws}/${site_id}"
    # vasp-inputs.sh names the Slurm job after the CIF basename (no extension),
    # and the self-resubmit chain keeps that name; track it for the wait loop.
    job_cif="${surface_id}_${adsorbate}_${site_id}.cif"
    job_name="${job_cif%.cif}"
    site_dirs+=("$dir")
    site_jobs+=("$job_name")

    # Skip re-submission when a terminal marker already exists (resume support).
    if [[ "$force" -eq 0 ]]; then
        for m in CONVERGED NOT_CONVERGED_MAX_GEN CRASHED SCF_ILL CHECK_ERROR; do
            if [[ -e "${dir}/${m}" ]]; then
                log "step 3: ${site_id} already finished (${m}); not resubmitting."
                continue 2
            fi
        done
    fi

    log "step 3: submitting VASP relaxation for ${site_id}"
    mkdir -p "$dir"
    cp -- "$init_cif" "${dir}/${job_cif}"
    cp -- "${script_dir}/check-convergence.sh" "${dir}/check-convergence.sh"
    chmod +x "${dir}/check-convergence.sh"
    (
        cd "$dir"
        bash "${script_dir}/vasp-inputs.sh" "$job_cif"
        echo 1 > .gen_count
        rm -f CONVERGED NOT_CONVERGED_MAX_GEN CRASHED SCF_ILL CHECK_ERROR
        sbatch --export=ALL,MAX_GEN="$max_gen" vasp-gam.slurm
    )
    submitted=$((submitted + 1))
done
log "step 3: ${submitted} job(s) submitted; waiting for completion markers ..."

# --- wait until every site directory carries a terminal marker ---------------
# A site is "settled" when it has a terminal marker. As a fallback for chains
# that die without writing one (hard node failure, wallclock kill), a site is
# also considered settled after its Slurm job has been absent from the queue
# for two consecutive polls with still no marker; Step 4 then skips it as
# unconverged. Two polls avoid a false positive during the brief gap while one
# generation exits and the next is submitted.
markers=(CONVERGED NOT_CONVERGED_MAX_GEN CRASHED SCF_ILL CHECK_ERROR)
declare -a dead_polls
for _ in "${site_dirs[@]}"; do dead_polls+=(0); done
while true; do
    pending=0
    for i in "${!site_dirs[@]}"; do
        dir=${site_dirs[$i]}
        done_one=0
        for m in "${markers[@]}"; do
            [[ -e "${dir}/${m}" ]] && { done_one=1; break; }
        done
        if [[ "$done_one" -eq 1 ]]; then
            dead_polls[$i]=0
            continue
        fi
        if job_alive "${site_jobs[$i]}"; then
            dead_polls[$i]=0
        else
            dead_polls[$i]=$(( ${dead_polls[$i]} + 1 ))
            if [[ "${dead_polls[$i]}" -ge 2 ]]; then
                log "step 3: WARNING: ${site_jobs[$i]} left the queue with no marker; treating as failed."
                continue
            fi
        fi
        pending=$((pending + 1))
    done
    [[ "$pending" -eq 0 ]] && break
    sleep "$poll"
done
log "step 3: all site relaxations settled."

# ============================================================================
# Step 4: compute adsorption energies and record them
# ============================================================================
recorded=0
skipped=0
for s in "${sites[@]}"; do
    site_id=$(site_id_of "$s")
    dir="${ads_ws}/${site_id}"

    marker=""
    for m in "${markers[@]}"; do
        [[ -e "${dir}/${m}" ]] && { marker=$m; break; }
    done
    if [[ "$marker" != "CONVERGED" ]]; then
        log "step 4: ${site_id} not converged (marker: ${marker:-none}); skipping."
        skipped=$((skipped + 1))
        continue
    fi

    contcar="${dir}/CONTCAR"
    outcar="${dir}/OUTCAR"
    [[ -s "$contcar" ]] || { log "step 4: ${site_id} missing CONTCAR; skipping."; skipped=$((skipped + 1)); continue; }
    [[ -s "$outcar" ]]  || { log "step 4: ${site_id} missing OUTCAR; skipping.";  skipped=$((skipped + 1)); continue; }

    relaxed_cif="${dir}/01_relaxed_adsorbate.cif"
    log "step 4: ${site_id} converting CONTCAR -> ${relaxed_cif}"
    "$PYTHON" - "$contcar" "$relaxed_cif" <<'PY'
import sys
from ase.io import read, write
write(sys.argv[2], read(sys.argv[1], format="vasp"))
PY

    # Final total energy (energy without entropy, matching workflow 1), then
    # E_ads = E(slab+ads) - E(slab) - E_ref.
    e_ads=$(
        "$PYTHON" - "$outcar" "$slab_energy" "$ref_energy" <<'PY'
import re, sys
total = None
with open(sys.argv[1]) as fh:
    for line in fh:
        m = re.search(r'energy\s+without\s+entropy\s*=\s*([+-]?\d+(?:\.\d+)?)', line)
        if m:
            total = float(m.group(1))
if total is None:
    print("")
else:
    print(repr(total - float(sys.argv[2]) - float(sys.argv[3])))
PY
    )
    if [[ -z "$e_ads" ]]; then
        log "step 4: ${site_id} could not extract total energy from OUTCAR; skipping."
        skipped=$((skipped + 1))
        continue
    fi
    log "step 4: ${site_id} E_ads = ${e_ads} eV"

    # A relaxed adsorbate that has moved beyond record-energy's validation
    # limit no longer aborts the whole multi-site workflow. Preserve all other
    # record-energy failures as fatal so filesystem/database problems are not
    # silently hidden.
    record_output=""
    if record_output=$("$PYTHON" "${repo_dir}/hea_dataset.py" record-energy \
        --root "$root" --surface-id "$surface_id" \
        --adsorbate "$adsorbate" --site-id "$s" \
        --relaxed-cif "$relaxed_cif" --energy "$e_ads" \
        --status computed \
        --notes "auto: E_slab=${slab_energy} eV, E_ref=${ref_energy} eV" 2>&1
    ); then
        [[ -z "$record_output" ]] || printf '%s\n' "$record_output" >&2
    else
        record_status=$?
        [[ -z "$record_output" ]] || printf '%s\n' "$record_output" >&2
        if [[ "$record_output" == *"Adsorbate moved too far during relaxation:"* ]]; then
            log "step 4: ${site_id} adsorbate moved too far; skipping."
            skipped=$((skipped + 1))
            continue
        fi
        exit "$record_status"
    fi
    recorded=$((recorded + 1))
done

log "done. surface-id=${surface_id} adsorbate=${adsorbate} recorded=${recorded} skipped=${skipped}"
