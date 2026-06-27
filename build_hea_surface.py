#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Build a High-Entropy Alloy (HEA) surface model.

Workflow:
    1. Use an FCC(111) facet as the template; take 4 atoms along each of the
       x, y and z directions to build the slab.
    2. Add a vacuum layer along z (default total thickness 15 Angstrom).
    3. From the input element types and ratios (randomly generated when not
       given), convert the ratios into integer atom counts and randomly
       substitute them onto the atomic sites to obtain an initial structure.
    4. Optimize the element distribution with the SQS (Special Quasirandom
       Structure) method so that, in terms of the chosen cluster correlation
       functions, it best approximates an ideal random solid solution.
    5. Save the final structure in CIF format.

Dependencies:
    - ase   (pip install ase)
    - icet  (pip install icet)        # used for SQS optimization
    - numpy

Examples:
    # Specify elements and ratios
    python build_hea_surface.py --elements Fe Co Ni Cr Mn --ratios 1 1 1 1 1

    # Only give elements; ratios are generated randomly
    python build_hea_surface.py --elements Fe Co Ni Cr Mn

    # Customize vacuum, supercell size, SQS steps and output file
    python build_hea_surface.py -e Cu Ni Pd -s 4 4 4 --vacuum 15 \
        --n-steps 20000 -o CuNiPd_111_sqs.cif
"""

import argparse
import sys

import numpy as np
from ase.build import fcc111
from ase.io import write

# ---------------------------------------------------------------------------
# Lattice constants (Angstrom) of common FCC metals, used to estimate the
# template lattice constant via Vegard's law. Elements not listed here fall
# back to DEFAULT_A.
# ---------------------------------------------------------------------------
FCC_LATTICE_CONSTANTS = {
    "Al": 4.05, "Ni": 3.52, "Cu": 3.61, "Pd": 3.89, "Ag": 4.09,
    "Pt": 3.92, "Au": 4.08, "Pb": 4.95, "Rh": 3.80, "Ir": 3.84,
    "Ga": 4.51, "In": 4.59, "Sn": 4.89,
    # The following are not FCC in their ground state; the values are common
    # equivalent-FCC lattice-constant estimates used as an approximate
    # template inside an HEA solid solution.
    "Fe": 3.59, "Co": 3.54, "Cr": 3.62, "Mn": 3.78, "Ti": 4.08,
    "V": 3.78, "Mo": 3.86, "W": 3.93, "Zn": 3.94,
}
DEFAULT_A = 3.6  # default lattice constant (Angstrom) when no element is known


def parse_args(argv=None):
    p = argparse.ArgumentParser(
        description="Build an FCC(111) HEA surface model and optimize the "
                    "element distribution with SQS.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    p.add_argument(
        "-e", "--elements", nargs="+", required=True, metavar="EL",
        help="Element types to include, e.g. Fe Co Ni Cr Mn",
    )
    p.add_argument(
        "-r", "--ratios", nargs="+", type=float, default=None, metavar="R",
        help="Ratio of each element (one per element in --elements). "
             "Randomly generated when not specified.",
    )
    p.add_argument(
        "-s", "--size", nargs=3, type=int, default=(4, 4, 4),
        metavar=("NX", "NY", "NZ"),
        help="Number of repeats/layers along x, y and z.",
    )
    p.add_argument(
        "--vacuum", type=float, default=15.0,
        help="Total vacuum-layer thickness (Angstrom), split symmetrically "
             "on both sides of the slab.",
    )
    p.add_argument(
        "-a", "--lattice-constant", type=float, default=None,
        help="Template lattice constant (Angstrom). Defaults to a Vegard's-law "
             "estimate from the elements.",
    )
    p.add_argument(
        "--cutoffs", nargs="+", type=float, default=[6.0, 4.5],
        help="Cluster-space cutoff radii (Angstrom) for icet: pair, triplet, ...",
    )
    p.add_argument(
        "--n-steps", type=int, default=10000,
        help="Number of Monte Carlo steps for SQS optimization.",
    )
    p.add_argument(
        "--seed", type=int, default=None,
        help="Random seed (for ratio generation, random substitution and SQS) "
             "for reproducibility.",
    )
    p.add_argument(
        "--no-sqs", action="store_true",
        help="Skip SQS optimization and output only the randomly substituted "
             "structure.",
    )
    p.add_argument(
        "-o", "--output", default=None,
        help="Output CIF filename. Auto-named from the elements by default.",
    )
    return p.parse_args(argv)


def normalize_ratios(elements, ratios, rng):
    """Return normalized ratios (summing to 1); generate them when not given."""
    n = len(elements)
    if ratios is None:
        # Sample uniformly on (0,1) then normalize -> random, non-zero ratios.
        ratios = rng.random(n) + 0.1
        print("[info] No ratios given; generated randomly.")
    else:
        if len(ratios) != n:
            sys.exit(f"[error] Number of ratios ({len(ratios)}) does not match "
                     f"number of elements ({n}).")
        ratios = np.asarray(ratios, dtype=float)
        if np.any(ratios < 0):
            sys.exit("[error] Ratios must not be negative.")
    total = ratios.sum()
    if total <= 0:
        sys.exit("[error] The sum of ratios must be greater than 0.")
    return ratios / total


def largest_remainder_counts(fractions, n_sites):
    """Allocate integer counts summing exactly to n_sites (largest remainder)."""
    raw = fractions * n_sites
    floor = np.floor(raw).astype(int)
    remainder = n_sites - floor.sum()
    if remainder > 0:
        # Give +1 to the elements with the largest fractional remainders.
        order = np.argsort(-(raw - floor))
        floor[order[:remainder]] += 1
    return floor


def estimate_lattice_constant(elements, fractions):
    """Estimate the template lattice constant via Vegard's law (ratio-weighted)."""
    a_vals = np.array(
        [FCC_LATTICE_CONSTANTS.get(el, DEFAULT_A) for el in elements]
    )
    return float(np.dot(fractions, a_vals))


def build_template_slab(symbol, size, a, vacuum_total):
    """Build the FCC(111) slab template.

    ASE's center(vacuum=v) leaves v of vacuum on each side of the slab, so the
    total vacuum thickness is 2v; passing vacuum_total / 2 makes the total
    thickness equal to the requested value.
    """
    slab = fcc111(symbol=symbol, size=size, a=a, vacuum=vacuum_total / 2.0)
    # SQS / icet need a 3D-periodic supercell; the vacuum keeps the z-direction
    # images isolated from each other.
    slab.set_pbc(True)
    return slab


def assign_random_substitution(slab, elements, counts, rng):
    """Randomly substitute the given count of each element onto the slab sites."""
    symbols = []
    for el, c in zip(elements, counts):
        symbols.extend([el] * int(c))
    rng.shuffle(symbols)
    slab = slab.copy()
    slab.set_chemical_symbols(symbols)
    return slab


def run_sqs(slab, elements, counts, cutoffs, n_steps, seed):
    """Generate an SQS structure with icet on the fixed supercell (slab)."""
    try:
        from icet import ClusterSpace
        from icet.tools.structure_generation import (
            generate_sqs_from_supercells,
        )
    except ImportError:
        print(
            "[warn] icet not installed; skipping SQS and outputting the random "
            "structure.\n"
            "       For SQS run: pip install icet",
            file=sys.stderr,
        )
        return None

    n_sites = len(slab)
    # Derive concentrations that are exactly realizable on this supercell from
    # the integer counts, avoiding non-divisible-concentration errors.
    target_conc = {el: int(c) / n_sites for el, c in zip(elements, counts)}

    # Use the slab as the parent structure of the cluster space, allowing every
    # element to occupy every site.
    cs = ClusterSpace(
        structure=slab,
        cutoffs=list(cutoffs),
        chemical_symbols=list(elements),
    )
    print(f"[info] Starting SQS optimization: {n_sites} sites, {n_steps} steps ...")
    sqs = generate_sqs_from_supercells(
        cluster_space=cs,
        supercells=[slab],
        target_concentrations=target_conc,
        n_steps=n_steps,
        random_seed=seed,
    )
    return sqs


def main(argv=None):
    args = parse_args(argv)
    rng = np.random.default_rng(args.seed)

    elements = args.elements
    fractions = normalize_ratios(elements, args.ratios, rng)

    size = tuple(args.size)
    a = args.lattice_constant or estimate_lattice_constant(elements, fractions)

    # 1) Build the template slab (first element as placeholder; substituted later).
    slab = build_template_slab(elements[0], size, a, args.vacuum)
    n_sites = len(slab)

    # 2) ratios -> integer counts (summing exactly to the number of sites).
    counts = largest_remainder_counts(fractions, n_sites)

    print("=" * 60)
    print(f"Template      : FCC(111), size={size}, a={a:.3f} A")
    print(f"Total sites   : {n_sites}")
    print(f"Vacuum layer  : {args.vacuum} A (total)")
    print("Composition   :")
    for el, frac, c in zip(elements, fractions, counts):
        print(f"  {el:>3s} : target {frac*100:6.2f}%  ->  {int(c):>3d} atoms "
              f"({int(c)/n_sites*100:6.2f}%)")
    print("=" * 60)

    # 3) Random substitution -> initial structure.
    random_slab = assign_random_substitution(slab, elements, counts, rng)

    # 4) SQS optimization (unless --no-sqs).
    final = None
    if not args.no_sqs:
        final = run_sqs(slab, elements, counts, args.cutoffs,
                        args.n_steps, args.seed)
    if final is None:
        final = random_slab
        method = "random"
    else:
        method = "sqs"
        print("[info] SQS optimization finished.")

    # 5) Save as CIF.
    out = args.output or f"HEA_{''.join(elements)}_111_{method}.cif"
    write(out, final)
    print(f"[done] Structure saved to: {out}  (method: {method})")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
