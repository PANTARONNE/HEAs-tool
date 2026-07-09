#!/usr/bin/env bash
#
# check-convergence.sh — inspect a finished VASP relaxation and report its state
# via an exit code. It is PURE: it only reads OUTCAR/OSZICAR/INCAR and never
# touches CONTCAR, never resubmits, and never deletes anything. All side effects
# (CONTCAR->POSCAR, resubmission, archiving) belong to the caller.
#
# Usage:
#   check-convergence.sh [calc_dir]      # defaults to the current directory
#
# Environment:
#   SCF_STALL_LIMIT   consecutive NELM-capped ionic steps that mark the SCF as
#                     ill-behaved (default 3).
#
# Exit codes:
#   0   converged      — VASP printed "reached required accuracy"
#   10  continue       — finished normally but not converged (e.g. hit NSW)
#   20  crashed        — no normal termination footer (timeout / segfault / OOM)
#   30  scf_ill        — SCF hit NELM for SCF_STALL_LIMIT consecutive ionic steps
#   2   usage/IO error

set -Eeuo pipefail

calc_dir=${1:-.}
outcar="${calc_dir}/OUTCAR"
oszicar="${calc_dir}/OSZICAR"
incar="${calc_dir}/INCAR"
scf_stall_limit=${SCF_STALL_LIMIT:-10}

incar_int() {
    # Read an integer INCAR tag (first match), fall back to a default.
    local tag=$1 default=$2 value
    value=$(
        awk -v tag="$tag" '
            BEGIN { IGNORECASE = 1 }
            {
                line = $0
                sub(/[#!].*/, "", line)          # strip inline comments
                if (match(line, "(^|[[:space:]])" tag "[[:space:]]*=[[:space:]]*[0-9]+")) {
                    split(line, parts, /=/)
                    gsub(/[^0-9]/, "", parts[2])
                    if (parts[2] != "") { print parts[2]; exit }
                }
            }
        ' "$incar" 2>/dev/null
    )
    [[ -n "$value" ]] && printf '%s\n' "$value" || printf '%s\n' "$default"
}

[[ -f "$outcar" ]] || {
    printf '[check] OUTCAR not found in %s -> crashed\n' "$calc_dir" >&2
    exit 20
}

# 1) Normal termination: VASP writes the timing footer only on a clean exit.
if ! grep -q 'General timing and accounting' "$outcar"; then
    printf '[check] no timing footer -> crashed/timeout\n' >&2
    exit 20
fi

# 2) Ionic convergence: VASP's own verdict against EDIFFG. Different VASP
#    builds spell the trailing text differently, so match the stable prefix.
if grep -qi 'reached required accuracy' "$outcar"; then
    printf '[check] reached required accuracy -> converged\n' >&2
    exit 0
fi

# 3) SCF health: scan every ionic step's electronic-iteration count. A single
#    step brushing NELM is normal early in a relaxation; a *run* of consecutive
#    steps that all hit NELM means the SCF cannot converge and continuing would
#    only waste walltime. We report the longest such consecutive run.
nelm=$(incar_int NELM 60)
max_stall=$(
    awk -v nelm="$nelm" '
        # An ionic step closes on the "N F= ..." line; the electronic-step
        # counter (elec) accumulated since the previous close is that step size.
        /^[[:space:]]*[0-9]+[[:space:]]+F=/ {
            if (nelm > 0 && elec >= nelm) {
                run++
                if (run > max_run) max_run = run
            } else {
                run = 0
            }
            elec = 0
            next
        }
        # Electronic SCF lines look like "DAV:  12   ..." or "RMM:  7  ...".
        /^[[:space:]]*[A-Za-z]+[[:space:]]*:[[:space:]]*[0-9]+/ { elec++ }
        END { printf "%d\n", max_run + 0 }
    ' "$oszicar" 2>/dev/null || echo 0
)

if [[ "${max_stall:-0}" -ge "$scf_stall_limit" && "$scf_stall_limit" -gt 0 ]]; then
    printf '[check] %s consecutive ionic steps hit NELM=%s -> SCF ill-behaved\n' \
        "$max_stall" "$nelm" >&2
    exit 30
fi

# 4) Finished cleanly, SCF healthy, but not converged -> needs continuation.
ionic=$(grep -c '^[[:space:]]*[0-9]\+[[:space:]]\+F=' "$oszicar" 2>/dev/null || echo 0)
printf '[check] finished without convergence (%s ionic steps) -> continue\n' \
    "${ionic:-0}" >&2
exit 10
