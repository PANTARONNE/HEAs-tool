#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Generate an OpenMX static (single-point) SCF input (.dat) file from a CIF.

The generated input is configured for a non-self-consistent-geometry, single
SCF calculation that writes the Hamiltonian (H) and overlap (S) matrices to
``<System.Name>.scfout`` via the ``HS.fileout on`` keyword. That .scfout file
is the standard OpenMX product consumed downstream (Wannier90 interface,
NEGF/transport, tight-binding analysis, ...). The text H/S can be dumped from
it with OpenMX's ``read_scfout`` / ``analysis_example`` utilities.

Workflow:
    1. Read the CIF with ASE (lattice vectors + fractional coordinates).
    2. Map every element onto its OpenMX PAO basis, pseudopotential (VPS) and
       number of valence electrons (OpenMX 2019 "standard" recommendation).
    3. Split the valence electrons into spin up/down initial occupations,
       optionally seeded with an initial magnetic moment per element.
    4. Choose a Monkhorst-Pack k-grid from a target k-density (or use --kgrid).
    5. Write a ready-to-run .dat file.

Dependencies:
    - ase    (pip install ase)
    - numpy

IMPORTANT before running OpenMX:
    * Point --data-path at your OpenMX DFT_DATA directory (e.g. .../DFT_DATA19).
    * The PAO/VPS names below follow the 2019 standard recommendation; verify
      they match the files actually present in your DFT_DATA version.

Examples:
    python cif_to_openmx.py HEA_FeCoNiCrMn_111_sqs.cif \
        --data-path /home/me/openmx3.9/DFT_DATA19

    python cif_to_openmx.py struct.cif -o struct.dat \
        --kgrid 4 4 1 --energycutoff 300 --moments Fe=2.5 Co=1.6 Ni=0.6
"""

import argparse
import os
import sys

import numpy as np

# ---------------------------------------------------------------------------
# PAO basis / VPS pseudopotential / valence electrons, matched to the files
# actually present in ./DFT_DATA19 (PAO/*.pao, VPS/*.vps; valences taken from
# DFT_DATA19/vps_info.txt). Note that in this library Cr and Mn ship ONLY the
# regular variant (Cr6.0/Mn6.0, Cr_PBE19/Mn_PBE19) -- there is no "H" version,
# whereas Fe/Co/Ni/Cu/Zn use the hard ("H") variants. The s3p2d1 contraction
# is a light "standard" choice; bump the d/f count (e.g. s3p2d2) for accuracy.
# Add more elements as needed, then re-run.
# ---------------------------------------------------------------------------
SPECIES_DB = {
    # light / main-group (subset)
    "H":  ("H6.0-s2p1",      "H_PBE19",    1.0),
    "B":  ("B7.0-s2p2d1",    "B_PBE19",    3.0),
    "C":  ("C6.0-s2p2d1",    "C_PBE19",    4.0),
    "N":  ("N6.0-s2p2d1",    "N_PBE19",    5.0),
    "O":  ("O6.0-s2p2d1",    "O_PBE19",    6.0),
    "F":  ("F6.0-s2p2d1",    "F_PBE19",    7.0),
    "Al": ("Al7.0-s2p2d1",   "Al_PBE19",   3.0),
    "Si": ("Si7.0-s2p2d1",   "Si_PBE19",   4.0),
    "P":  ("P7.0-s2p2d1",    "P_PBE19",    5.0),
    "S":  ("S7.0-s2p2d1",    "S_PBE19",    6.0),
    # 3d transition metals (matched to DFT_DATA19)
    "Ti": ("Ti7.0-s3p2d1",   "Ti_PBE19",  12.0),
    "V":  ("V6.0-s3p2d1",    "V_PBE19",   13.0),
    "Cr": ("Cr6.0-s3p2d1",   "Cr_PBE19",  14.0),  # no H variant in this library
    "Mn": ("Mn6.0-s3p2d1",   "Mn_PBE19",  15.0),  # no H variant in this library
    "Fe": ("Fe6.0H-s3p2d1",  "Fe_PBE19H", 16.0),
    "Co": ("Co6.0H-s3p2d1",  "Co_PBE19H", 17.0),
    "Ni": ("Ni6.0H-s3p2d1",  "Ni_PBE19H", 18.0),
    "Cu": ("Cu6.0H-s3p2d1",  "Cu_PBE19H", 19.0),
    "Zn": ("Zn6.0H-s3p2d1",  "Zn_PBE19H", 20.0),
    # 4d/5d transition metals and post-transition metals requested for HEAs.
    "Mo": ("Mo7.0-s3p2d1",   "Mo_PBE19",  14.0),
    "W":  ("W7.0-s3p2d1",    "W_PBE19",   12.0),
    "Ga": ("Ga7.0-s2p2d1",   "Ga_PBE19",  13.0),
    "In": ("In7.0-s2p2d1",   "In_PBE19",  13.0),
    "Sn": ("Sn7.0-s2p2d1",   "Sn_PBE19",  14.0),
}

# Rough ferromagnetic initial moments (mu_B) used only to seed the SCF; the
# self-consistent loop relaxes them. Elements not listed start non-magnetic.
DEFAULT_MOMENTS = {
    "Cr": 1.0, "Mn": 3.0, "Fe": 3.0, "Co": 2.0, "Ni": 1.0,
}


def parse_args(argv=None):
    p = argparse.ArgumentParser(
        description="Generate an OpenMX static SCF input (.dat) for H/S matrix "
                    "output (HS.fileout on) from a CIF structure.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    p.add_argument("cif", help="Input structure in CIF format.")
    p.add_argument(
        "-o", "--output", default=None,
        help="Output .dat filename. Defaults to <cif-stem>.dat.",
    )
    p.add_argument(
        "--data-path", required=True,
        help="REQUIRED. Server-side OpenMX DFT_DATA directory (where OpenMX will "
             "run), written verbatim as DATA.PATH in the .dat file, "
             "e.g. /home/user/openmx3.9/DFT_DATA19.",
    )
    p.add_argument(
        "--xc", default="GGA-PBE",
        help="Exchange-correlation functional (scf.XcType), e.g. GGA-PBE, LDA.",
    )
    p.add_argument(
        "--no-spin", action="store_true",
        help="Disable spin polarization (scf.SpinPolarization off). The "
             "FeCoNiCrMn-type alloys are magnetic, so spin-on is the default.",
    )
    p.add_argument(
        "--kgrid", nargs=3, type=int, default=None, metavar=("KX", "KY", "KZ"),
        help="Explicit Monkhorst-Pack grid. Auto-chosen from --kdens otherwise.",
    )
    p.add_argument(
        "--kdens", type=float, default=30.0,
        help="Target k-density length R (Angstrom) for the auto k-grid: "
             "n_i = round(R / |a_i|), min 1. Larger R -> denser grid.",
    )
    p.add_argument(
        "--energycutoff", type=float, default=300.0,
        help="Real-space grid cutoff scf.energycutoff (Ry).",
    )
    p.add_argument(
        "--maxiter", type=int, default=200,
        help="Maximum number of SCF iterations (scf.maxIter).",
    )
    p.add_argument(
        "--criterion", type=float, default=1.0e-6,
        help="SCF convergence criterion in Hartree (scf.criterion).",
    )
    p.add_argument(
        "--eltemp", type=float, default=300.0,
        help="Electronic temperature in K (scf.ElectronicTemperature).",
    )
    p.add_argument(
        "--moments", nargs="*", default=None, metavar="EL=MU",
        help="Override initial magnetic moments, e.g. Fe=2.5 Co=1.6 Cr=-1.0. "
             "Used only with spin on.",
    )
    p.add_argument(
        "--coord-unit", choices=["FRAC", "Ang"], default="FRAC",
        help="Unit for Atoms.SpeciesAndCoordinates.",
    )
    return p.parse_args(argv)


def read_structure(path):
    """Read a CIF into an ASE Atoms object."""
    try:
        from ase.io import read
    except ImportError:
        sys.exit("[error] ASE is required: pip install ase")
    if not os.path.isfile(path):
        sys.exit(f"[error] CIF file not found: {path}")
    atoms = read(path)
    if atoms.cell.rank != 3:
        sys.exit("[error] The structure has no full 3D cell; OpenMX needs "
                 "Atoms.UnitVectors. Check the CIF.")
    return atoms


def parse_moment_overrides(items):
    """Parse ['Fe=2.5', 'Co=1.6'] into {'Fe': 2.5, 'Co': 1.6}."""
    moments = dict(DEFAULT_MOMENTS)
    if not items:
        return moments
    for item in items:
        if "=" not in item:
            sys.exit(f"[error] Bad --moments entry '{item}', expected EL=MU.")
        el, val = item.split("=", 1)
        try:
            moments[el] = float(val)
        except ValueError:
            sys.exit(f"[error] Bad moment value in '{item}'.")
    return moments


def auto_kgrid(cell, kdens):
    """Pick a Monkhorst-Pack grid from a target k-density length (Angstrom)."""
    lengths = np.linalg.norm(np.asarray(cell), axis=1)
    return [max(1, int(round(kdens / L))) for L in lengths]


def spin_occupation(valence, moment, spin_on):
    """Split valence electrons into (up, down) initial occupations."""
    if not spin_on:
        half = valence / 2.0
        return half, half
    m = max(-valence, min(valence, moment))  # clamp to physical range
    up = (valence + m) / 2.0
    down = (valence - m) / 2.0
    return up, down


def build_dat(atoms, args, moments):
    """Assemble the full OpenMX .dat content as a string."""
    spin_on = not args.no_spin
    symbols = atoms.get_chemical_symbols()
    species = sorted(set(symbols))

    # Validate that every element is known.
    missing = [s for s in species if s not in SPECIES_DB]
    if missing:
        sys.exit(f"[error] No PAO/VPS entry for: {', '.join(missing)}. "
                 f"Add them to SPECIES_DB.")

    system_name = os.path.splitext(os.path.basename(args.output))[0]
    cell = atoms.get_cell()[:]
    kgrid = args.kgrid or auto_kgrid(cell, args.kdens)

    if args.coord_unit == "FRAC":
        coords = atoms.get_scaled_positions(wrap=True)
    else:
        coords = atoms.get_positions()

    lines = []
    add = lines.append

    # --- header ------------------------------------------------------------
    add(f"# OpenMX static SCF input auto-generated from {os.path.basename(args.cif)}")
    add(f"# Formula : {atoms.get_chemical_formula()}")
    add(f"# Atoms   : {len(atoms)}   Species: {len(species)}")
    add("# Purpose : write Hamiltonian/overlap matrices (HS.fileout on).")
    add("")

    # --- files / verbosity -------------------------------------------------
    add("System.CurrentDirectory     ./")
    add(f"System.Name                 {system_name}")
    add(f"DATA.PATH                   {args.data_path}")
    add("level.of.stdout             1")
    add("level.of.fileout            1")
    add("")

    # --- species -----------------------------------------------------------
    add(f"Species.Number              {len(species)}")
    add("<Definition.of.Atomic.Species")
    for s in species:
        pao, vps, _ = SPECIES_DB[s]
        add(f"  {s:<3s} {pao:<16s} {vps}")
    add("Definition.of.Atomic.Species>")
    add("")

    # --- atoms / coordinates ----------------------------------------------
    add(f"Atoms.Number                {len(atoms)}")
    add(f"Atoms.SpeciesAndCoordinates.Unit   {args.coord_unit}")
    add("<Atoms.SpeciesAndCoordinates")
    for i, (s, xyz) in enumerate(zip(symbols, coords), start=1):
        _, _, val = SPECIES_DB[s]
        up, down = spin_occupation(val, moments.get(s, 0.0), spin_on)
        add(f"  {i:>4d} {s:<3s} "
            f"{xyz[0]:18.12f} {xyz[1]:18.12f} {xyz[2]:18.12f} "
            f"{up:8.4f} {down:8.4f}")
    add("Atoms.SpeciesAndCoordinates>")
    add("")

    add("Atoms.UnitVectors.Unit      Ang")
    add("<Atoms.UnitVectors")
    for v in cell:
        add(f"  {v[0]:18.12f} {v[1]:18.12f} {v[2]:18.12f}")
    add("Atoms.UnitVectors>")
    add("")

    # --- SCF parameters ----------------------------------------------------
    add(f"scf.XcType                  {args.xc}")
    add(f"scf.SpinPolarization        {'on' if spin_on else 'off'}")
    add(f"scf.ElectronicTemperature   {args.eltemp}")
    add(f"scf.energycutoff            {args.energycutoff}")
    add(f"scf.maxIter                 {args.maxiter}")
    add("scf.EigenvalueSolver        Band")
    add(f"scf.Kgrid                   {kgrid[0]} {kgrid[1]} {kgrid[2]}")
    add("scf.Mixing.Type             rmm-diisk")
    add("scf.Init.Mixing.Weight      0.30")
    add("scf.Min.Mixing.Weight       0.001")
    add("scf.Max.Mixing.Weight       0.400")
    add("scf.Mixing.History          7")
    add("scf.Mixing.StartPulay       6")
    add(f"scf.criterion               {args.criterion}")
    add("")

    # --- single point, no geometry optimization ---------------------------
    add("MD.Type                     Nomd")
    add("MD.maxIter                  1")
    add("")

    # --- Hamiltonian / overlap output -------------------------------------
    add("HS.fileout                  on")
    add("")

    return "\n".join(lines) + "\n", system_name, species, kgrid


def main(argv=None):
    args = parse_args(argv)
    if args.output is None:
        stem = os.path.splitext(os.path.basename(args.cif))[0]
        # System.Name must be a plain token; keep it filesystem-friendly.
        args.output = "".join(c if c.isalnum() or c in "._-" else "_"
                              for c in stem) + ".dat"

    moments = parse_moment_overrides(args.moments)
    atoms = read_structure(args.cif)
    content, system_name, species, kgrid = build_dat(atoms, args, moments)

    with open(args.output, "w") as f:
        f.write(content)

    print("=" * 64)
    print(f"Input written : {args.output}")
    print(f"System.Name   : {system_name}")
    print(f"Atoms         : {len(atoms)}   Species: {', '.join(species)}")
    print(f"Spin          : {'off' if args.no_spin else 'on'}")
    print(f"k-grid        : {kgrid[0]} {kgrid[1]} {kgrid[2]}")
    print(f"Energy cutoff : {args.energycutoff} Ry")
    print("=" * 64)
    print("Run with, e.g.:  mpirun -np <N> openmx "
          f"{args.output} > {system_name}.std")
    print(f"H/S matrices will be written to: {system_name}.scfout")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
