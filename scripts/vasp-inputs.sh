#!/usr/bin/env bash

set -Eeuo pipefail

usage() {
    printf 'Usage: %s <structure.cif>\n' "$(basename "$0")" >&2
    exit 2
}

[[ $# -eq 1 ]] || usage

cif_file=$1
script_dir=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)
incar_template="${script_dir}/INCAR-opt"
slurm_template="${script_dir}/vasp-gam.slurm"

[[ -f "$cif_file" ]] || {
    printf 'Error: CIF file not found: %s\n' "$cif_file" >&2
    exit 1
}
[[ -f "$incar_template" ]] || {
    printf 'Error: INCAR template not found: %s\n' "$incar_template" >&2
    exit 1
}
[[ -f "$slurm_template" ]] || {
    printf 'Error: Slurm template not found: %s\n' "$slurm_template" >&2
    exit 1
}
command -v vaspkit >/dev/null 2>&1 || {
    printf 'Error: vaspkit is not available in PATH.\n' >&2
    exit 1
}

# Read the elements from the CIF atom-site loop, preserving first occurrence.
element_order=$(
    awk '
        BEGIN { in_loop = 0; field_count = 0; symbol_field = 0 }

        /^[[:space:]]*loop_[[:space:]]*$/ {
            in_loop = 1
            field_count = 0
            symbol_field = 0
            next
        }

        in_loop && /^[[:space:]]*_/ {
            field_count++
            field_name = $1
            sub(/\r$/, "", field_name)
            if (field_name == "_atom_site_type_symbol") {
                symbol_field = field_count
            }
            next
        }

        in_loop && symbol_field && /^[[:space:]]*#/ { next }
        in_loop && symbol_field && /^[[:space:]]*$/ { next }

        in_loop && symbol_field {
            symbol = $symbol_field
            gsub(/^["\047]|["\047\r]$/, "", symbol)
            if (symbol != "" && symbol != "." && symbol != "?" && !seen[symbol]++) {
                order[++count] = symbol
            }
        }

        END {
            for (i = 1; i <= count; i++) {
                printf "%s%s", (i == 1 ? "" : " "), order[i]
            }
        }
    ' "$cif_file"
)

[[ -n "$element_order" ]] || {
    printf 'Error: could not read _atom_site_type_symbol from CIF: %s\n' "$cif_file" >&2
    exit 1
}

# 105: explicitly pass the CIF element order instead of relying on its default.
printf '105\n%s\n%s\n' "$cif_file" "$element_order" | vaspkit > /dev/null 2>&1
[[ -s POSCAR ]] || {
    printf 'Error: VASPKIT task 105 did not generate POSCAR.\n' >&2
    exit 1
}

poscar_order=$(awk 'NR == 6 { for (i = 1; i <= NF; i++) printf "%s%s", (i == 1 ? "" : " "), $i }' POSCAR)
[[ "$poscar_order" == "$element_order" ]] || {
    printf 'Error: POSCAR element order differs from CIF.\n' >&2
    printf '  CIF:    %s\n  POSCAR: %s\n' "$element_order" "$poscar_order" >&2
    exit 1
}

cp -- "$incar_template" INCAR

# 103: default POTCAR. 102: M-P scheme, Gamma-only.
printf '103\n' | vaspkit > /dev/null 2>&1
printf '102\n1\n0\n' | vaspkit > /dev/null 2>&1

# 402: POSCAR -> fix atoms by z range -> fractional z in [0, 0.46].
# The resulting F F F flags constrain motion in all three Cartesian directions.
printf '402\n1\n3\n0 0.46\n1\nall\n' | vaspkit > /dev/null 2>&1
[[ -s POSCAR_FIX.vasp ]] || {
    printf 'Error: VASPKIT task 402 did not generate POSCAR_FIX.vasp\n' >&2
    exit 1
}
mv -- POSCAR_FIX.vasp POSCAR

for output in INCAR POSCAR POTCAR KPOINTS; do
    [[ -s "$output" ]] || {
        printf 'Error: expected output is missing or empty: %s\n' "$output" >&2
        exit 1
    }
done

cif_name=$(basename -- "$cif_file")
job_name=${cif_name%.[cC][iI][fF]}
[[ "$job_name" != *[[:space:]]* ]] || {
    printf 'Error: the CIF basename cannot contain whitespace in a Slurm job name: %s\n' "$job_name" >&2
    exit 1
}

cp -- "$slurm_template" vasp-gam.slurm
awk -v job_name="$job_name" '
    {
        sub(/\r$/, "")
    }
    NR == 1 {
        print
        print "#SBATCH --job-name=" job_name
        next
    }
    /^[[:space:]]*#SBATCH[[:space:]]+(-J|--job-name)([=[:space:]]|$)/ {
        next
    }
    {
        print
    }
' vasp-gam.slurm > vasp-gam.slurm.tmp
mv -- vasp-gam.slurm.tmp vasp-gam.slurm

printf 'Generated VASP inputs and vasp-gam.slurm (job name: %s) in %s\n' "$job_name" "$PWD"
