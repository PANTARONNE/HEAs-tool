#!/usr/bin/env bash
#
# Compute and record the adsorption energy for one already-calculated site.
# Unlike run-adsorption-workflow.sh, this script does not require a CONVERGED
# marker, so it can be used after manually recovering a site marked CRASHED.
#
# Usage:
#   record-single-site-adsorption.sh SURFACE_ID SITE [options]
#
# Examples:
#   bash scripts/record-single-site-adsorption.sh \
#     Co_13-Cr_13-Fe_13-Mn_12-Ni_13 5
#   bash scripts/record-single-site-adsorption.sh SURFACE_ID 5 -a NH2
#
# Options:
#   -a, --adsorbate SP    adsorbate species (auto-detected when unambiguous)
#       --root DIR        dataset root (default: dataset)
#       --workspace DIR   calculation workspace (default: workspace)
#       --calc-dir DIR    site calculation directory (default: inferred)
#       --slab-energy E   override E(slab) in eV (default: dataset index)
#       --ref-energy E    override E(ref) in eV (default: table below)
#       --max-adsorbate-displacement A
#                        validation limit passed to record-energy (default: 2.0)
#       --status STATUS   recorded energy status (default: computed)
#       --notes TEXT      append text to the generated provenance note
#   -h, --help
#
# The expected default calculation directory is:
#   workspace/SURFACE_ID/adsorbate/SP/site_XXXX/{CONTCAR,OUTCAR}

set -Eeuo pipefail

declare -A ADSORBATE_REF_ENERGY=(
    [N]="-8.3177"
    [NH]="-8.1023"
    [NH2]="-13.5359"
    [NH3]="-19.5405"
)

script_dir=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)
repo_dir=$(cd -- "${script_dir}/.." && pwd)
PYTHON=${PYTHON:-python}

root="dataset"
workspace="workspace"
adsorbate=""
calc_dir=""
slab_energy=""
ref_energy=""
max_displacement="2.0"
status="computed"
extra_notes=""
positionals=()

die() { printf 'Error: %s\n' "$*" >&2; exit 1; }
log() { printf '[single-site-energy] %s\n' "$*" >&2; }
need_value() { [[ $# -ge 2 ]] || die "$1 requires a value."; }

while [[ $# -gt 0 ]]; do
    case "$1" in
        -a|--adsorbate)
            need_value "$@"; adsorbate=$2; shift 2 ;;
        --root)
            need_value "$@"; root=$2; shift 2 ;;
        --workspace)
            need_value "$@"; workspace=$2; shift 2 ;;
        --calc-dir)
            need_value "$@"; calc_dir=$2; shift 2 ;;
        --slab-energy)
            need_value "$@"; slab_energy=$2; shift 2 ;;
        --ref-energy)
            need_value "$@"; ref_energy=$2; shift 2 ;;
        --max-adsorbate-displacement)
            need_value "$@"; max_displacement=$2; shift 2 ;;
        --status)
            need_value "$@"; status=$2; shift 2 ;;
        --notes)
            need_value "$@"; extra_notes=$2; shift 2 ;;
        -h|--help)
            awk 'NR==1{next} /^#/{sub(/^# ?/,""); print; next} {exit}' "$0"
            exit 0
            ;;
        --)
            shift
            while [[ $# -gt 0 ]]; do positionals+=("$1"); shift; done
            ;;
        -*) die "unknown argument: $1" ;;
        *) positionals+=("$1"); shift ;;
    esac
done

[[ ${#positionals[@]} -eq 2 ]] || die \
    "usage: $0 SURFACE_ID SITE [options]"
surface_id=${positionals[0]}
site=${positionals[1]}
[[ "$site" =~ ^[0-9]+$ ]] || die "SITE must be a positive integer."
site=$((10#$site))
[[ "$site" -ge 1 ]] || die "SITE must be a positive integer."
printf -v site_id 'site_%04d' "$site"

[[ -f "${repo_dir}/hea_dataset.py" ]] || die \
    "hea_dataset.py not found in ${repo_dir}."
command -v "$PYTHON" >/dev/null 2>&1 || die "Python command not found: ${PYTHON}"

# With the normal directory layout, SURFACE_ID and SITE are sufficient: infer
# the species from calculation/data directories. Require -a only if more than
# one adsorbate species exists for the same site.
if [[ -z "$adsorbate" ]]; then
    declare -A species_seen=()
    shopt -s nullglob
    for candidate in \
        "${workspace}/${surface_id}/adsorbate"/*/"${site_id}" \
        "${root}/${surface_id}/adsorbates"/*/"${site_id}"; do
        [[ -d "$candidate" ]] || continue
        species_dir=$(dirname -- "$candidate")
        species_seen["$(basename -- "$species_dir")"]=1
    done
    if [[ -n "$calc_dir" && "$(basename -- "$calc_dir")" == "$site_id" ]]; then
        species_seen["$(basename -- "$(dirname -- "$calc_dir")")"]=1
    fi
    shopt -u nullglob
    species=()
    for detected_species in "${!species_seen[@]}"; do
        species+=("$detected_species")
    done
    if [[ ${#species[@]} -eq 1 ]]; then
        adsorbate=${species[0]}
        log "auto-detected adsorbate: ${adsorbate}"
    elif [[ ${#species[@]} -eq 0 ]]; then
        die "could not infer the adsorbate for ${surface_id}/${site_id}; pass -a SP."
    else
        die "multiple adsorbates found for ${surface_id}/${site_id} (${species[*]}); pass -a SP."
    fi
fi
case "$adsorbate" in
    N|NH|NH2|NH3) ;;
    *) die "unsupported adsorbate: ${adsorbate} (choose N, NH, NH2 or NH3)." ;;
esac

if [[ -z "$calc_dir" ]]; then
    calc_dir="${workspace}/${surface_id}/adsorbate/${adsorbate}/${site_id}"
fi
contcar="${calc_dir}/CONTCAR"
outcar="${calc_dir}/OUTCAR"
[[ -s "$contcar" ]] || die "missing or empty CONTCAR: ${contcar}"
[[ -s "$outcar" ]] || die "missing or empty OUTCAR: ${outcar}"

for marker in CRASHED SCF_ILL CHECK_ERROR NOT_CONVERGED_MAX_GEN; do
    if [[ -e "${calc_dir}/${marker}" ]]; then
        log "WARNING: found ${marker}; continuing because this script intentionally ignores workflow markers."
        break
    fi
done

energy_record="${root}/${surface_id}/adsorbates/${adsorbate}/${site_id}/adsorption_energy.json"
[[ -f "$energy_record" ]] || die \
    "adsorption record not found: ${energy_record}. Run create-adsorbate-records first."

if [[ -z "$ref_energy" ]]; then
    ref_energy=${ADSORBATE_REF_ENERGY[$adsorbate]:-}
fi
[[ -n "$ref_energy" ]] || die \
    "no reference energy for ${adsorbate}; pass --ref-energy."

if [[ -z "$slab_energy" ]]; then
    log "reading E(slab) from the dataset index ..."
    slab_energy=$(
        "$PYTHON" - "$root" "$surface_id" <<'PY'
import sqlite3
import sys
from pathlib import Path

root, surface_id = sys.argv[1], sys.argv[2]
db = Path(root) / "index.sqlite"
value = ""
if db.is_file():
    with sqlite3.connect(db) as con:
        columns = {
            row[1] for row in con.execute("PRAGMA table_info(surfaces)")
        }
        if "total_energy_eV" in columns:
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
slab_energy=${slab_energy%$'\r'}
[[ -n "$slab_energy" ]] || die \
    "E(slab) is not recorded for ${surface_id}; pass --slab-energy."

readarray -t energies < <(
    "$PYTHON" - "$outcar" "$slab_energy" "$ref_energy" <<'PY'
import re
import sys

number = r"[+-]?(?:\d+(?:\.\d*)?|\.\d+)(?:[Ee][+-]?\d+)?"
total = None
with open(sys.argv[1], encoding="utf-8", errors="replace") as fh:
    for line in fh:
        match = re.search(rf"energy\s+without\s+entropy\s*=\s*({number})", line)
        if match:
            total = float(match.group(1))
if total is None:
    raise SystemExit("No 'energy without entropy' value found in OUTCAR.")
slab = float(sys.argv[2])
reference = float(sys.argv[3])
print(repr(total))
print(repr(total - slab - reference))
PY
)
[[ ${#energies[@]} -eq 2 ]] || die "could not extract the final energy from ${outcar}."
total_energy=${energies[0]%$'\r'}
adsorption_energy=${energies[1]%$'\r'}

relaxed_cif="${calc_dir}/01_relaxed_adsorbate.cif"
log "converting CONTCAR -> ${relaxed_cif}"
"$PYTHON" - "$contcar" "$relaxed_cif" <<'PY'
import sys
from ase.io import read, write

write(sys.argv[2], read(sys.argv[1], format="vasp"))
PY

notes="manual single-site recovery: E_total=${total_energy} eV, E_slab=${slab_energy} eV, E_ref=${ref_energy} eV"
if [[ -n "$extra_notes" ]]; then
    notes="${notes}; ${extra_notes}"
fi

log "${adsorbate}/${site_id}: E_ads = ${adsorption_energy} eV"
"$PYTHON" "${repo_dir}/hea_dataset.py" record-energy \
    --root "$root" --surface-id "$surface_id" \
    --adsorbate "$adsorbate" --site-id "$site" \
    --relaxed-cif "$relaxed_cif" --energy "$adsorption_energy" \
    --max-adsorbate-displacement "$max_displacement" \
    --status "$status" --notes "$notes"

log "recorded ${surface_id}/${adsorbate}/${site_id}."
