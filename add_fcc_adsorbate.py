#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Add N/NHx adsorbates on FCC hollow sites of an HEA FCC(111) slab CIF.

The script detects the top or bottom surface layers from the slab normal. For
an FCC(111) ABC-stacked slab, FCC hollow sites above the top layer are tagged by
the in-plane positions of atoms in the third layer below the surface. By
default, the actual adsorbate position is then relaxed to the local hollow
center formed by the three nearest top-layer atoms, which is more robust for
relaxed slabs with rumpling or lateral atom displacements.

Examples:
    python add_fcc_adsorbate.py FeCoNiCrMn11111-org.cif --adsorbate N
    python add_fcc_adsorbate.py FeCoNiCrMn11111-org.cif --adsorbate NH2 --site 5
    python add_fcc_adsorbate.py FeCoNiCrMn11111-org.cif --adsorbate NH3 --site all
"""

import argparse
import math
import os
import sys

import numpy as np


METALS_DEFAULT = {"Fe", "Co", "Ni", "Cr", "Mn", "Cu", "Pd", "Pt", "Rh", "Ir"}
ADSORBATES = {"N", "NH", "NH2", "NH3"}


def parse_args(argv=None):
    p = argparse.ArgumentParser(
        description="Read an HEA slab CIF and add N/NHx on FCC hollow sites.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    p.add_argument("cif", help="Input HEA slab CIF.")
    p.add_argument(
        "-a", "--adsorbate", choices=sorted(ADSORBATES), default="N",
        help="Adsorbate to place. The N atom is the atom bound to the surface.",
    )
    p.add_argument(
        "-s", "--site", default="1",
        help="1-based FCC site index to use, or 'all' to write one CIF per site.",
    )
    p.add_argument(
        "--side", choices=["top", "bottom"], default="top",
        help="Surface side on which to place the adsorbate.",
    )
    p.add_argument(
        "--height", type=float, default=1.25,
        help="N height above the FCC hollow plane along the surface normal (Ang).",
    )
    p.add_argument(
        "--nh", type=float, default=1.02,
        help="N-H bond length for NHx fragments (Ang).",
    )
    p.add_argument(
        "--layer-tol", type=float, default=0.60,
        help="Tolerance for grouping slab atoms into layers along the normal (Ang).",
    )
    p.add_argument(
        "--site-mode", choices=["relaxed", "projected"], default="relaxed",
        help="'relaxed' uses third-layer atoms only to tag FCC sites, then places "
             "the site at the centroid of the three nearest top-layer atoms. "
             "'projected' uses the ideal third-layer projection directly.",
    )
    p.add_argument(
        "--metals", nargs="+", default=sorted(METALS_DEFAULT),
        help="Elements treated as slab atoms when detecting surface layers.",
    )
    p.add_argument(
        "-o", "--output", default=None,
        help="Output CIF path. Defaults to <input>_<adsorbate>_fcc<site>.cif. "
             "Ignored for --site all, where a numbered set is written.",
    )
    p.add_argument(
        "--list-sites", action="store_true",
        help="Only print detected FCC hollow sites and do not write a CIF.",
    )
    return p.parse_args(argv)


def read_cif(path):
    try:
        from ase.io import read
    except ImportError:
        sys.exit("[error] ASE is required: pip install ase")
    if not os.path.isfile(path):
        sys.exit(f"[error] CIF file not found: {path}")
    atoms = read(path)
    if atoms.cell.rank != 3:
        sys.exit("[error] The input structure must have a full 3D cell.")
    return atoms


def write_cif(path, atoms):
    try:
        from ase.io import write
    except ImportError:
        sys.exit("[error] ASE is required: pip install ase")
    atoms_to_write = atoms.copy()
    atoms_to_write.info.pop("occupancy", None)
    write(path, atoms_to_write)


def normal_from_cell(cell, side):
    """Return the outward slab normal from a cell whose c vector is slab-normal."""
    normal = np.asarray(cell[2], dtype=float)
    norm = np.linalg.norm(normal)
    if norm < 1.0e-12:
        sys.exit("[error] Cell c vector is zero; cannot define slab normal.")
    normal = normal / norm
    if side == "bottom":
        normal = -normal
    return normal


def group_layers(atoms, normal, metal_symbols, tol):
    symbols = np.array(atoms.get_chemical_symbols())
    positions = atoms.get_positions()
    metal_mask = np.array([sym in metal_symbols for sym in symbols])
    if not np.any(metal_mask):
        sys.exit("[error] No slab atoms found. Check --metals.")

    metal_indices = np.where(metal_mask)[0]
    heights = positions[metal_indices] @ normal
    order = np.argsort(heights)

    layers = []
    for idx in metal_indices[order]:
        h = float(positions[idx] @ normal)
        if not layers or abs(h - layers[-1]["mean"]) > tol:
            layers.append({"indices": [int(idx)], "mean": h})
        else:
            layers[-1]["indices"].append(int(idx))
            vals = [float(positions[i] @ normal) for i in layers[-1]["indices"]]
            layers[-1]["mean"] = float(np.mean(vals))

    if len(layers) < 3:
        sys.exit(
            f"[error] Need at least 3 slab layers for FCC sites; found {len(layers)}."
        )
    return layers


def canonical_fractional_xy(frac_xy, ndigits=10):
    wrapped = np.mod(frac_xy[:2], 1.0)
    wrapped[np.isclose(wrapped, 1.0, atol=10.0 ** -ndigits)] = 0.0
    return tuple(np.round(wrapped, ndigits))


def inplane_vector_from_frac_delta(delta_xy, cell, normal):
    vec = float(delta_xy[0]) * cell[0] + float(delta_xy[1]) * cell[1]
    return vec - np.dot(vec, normal) * normal


def nearest_surface_triangle(atoms, surface_indices, anchor_xy, cell, normal):
    scaled = atoms.get_scaled_positions(wrap=True)
    distances = []
    for idx in surface_indices:
        delta = scaled[idx, :2] - anchor_xy
        delta -= np.round(delta)
        inplane = inplane_vector_from_frac_delta(delta, cell, normal)
        distances.append((float(np.linalg.norm(inplane)), int(idx), delta))

    distances.sort(key=lambda item: item[0])
    nearest = distances[:3]
    if len(nearest) < 3:
        sys.exit("[error] Need at least 3 top-layer atoms to define a hollow site.")

    top_atoms = [idx for _, idx, _ in nearest]
    unwrapped_xy = np.array([anchor_xy + delta for _, _, delta in nearest])
    center_xy = np.mod(np.mean(unwrapped_xy, axis=0), 1.0)
    local_height = float(np.mean([atoms.positions[idx] @ normal for idx in top_atoms]))
    return center_xy, local_height, top_atoms


def point_on_plane_from_frac_xy(frac_xy, cell, normal, plane_height):
    point = np.array([float(frac_xy[0]), float(frac_xy[1]), 0.0]) @ cell
    height_now = float(point @ normal)
    return point + (plane_height - height_now) * normal


def detect_fcc_sites(atoms, normal, side, metal_symbols, layer_tol, site_mode):
    layers = group_layers(atoms, normal, metal_symbols, layer_tol)
    surface = layers[-1]
    third_layer = layers[-3]
    surface_height = surface["mean"]

    cell = np.asarray(atoms.cell)
    sites = []
    seen = set()
    scaled = atoms.get_scaled_positions(wrap=True)

    for atom_index in third_layer["indices"]:
        frac = scaled[atom_index].copy()
        key = canonical_fractional_xy(frac)
        if key in seen:
            continue
        seen.add(key)

        anchor_xy = np.array([key[0], key[1]])
        if site_mode == "relaxed":
            site_xy, plane_height, top_atoms = nearest_surface_triangle(
                atoms, surface["indices"], anchor_xy, cell, normal
            )
        else:
            site_xy = anchor_xy
            plane_height = surface_height
            top_atoms = []

        base = point_on_plane_from_frac_xy(site_xy, cell, normal, plane_height)

        sites.append({
            "anchor_frac_xy": key,
            "frac_xy": canonical_fractional_xy(site_xy),
            "plane_pos": base,
            "source_atom": atom_index,
            "top_atoms": top_atoms,
        })

    sites.sort(key=lambda item: (item["anchor_frac_xy"][1], item["anchor_frac_xy"][0]))
    return sites, layers


def make_tangent_basis(normal):
    trial = np.array([1.0, 0.0, 0.0])
    if abs(float(np.dot(trial, normal))) > 0.90:
        trial = np.array([0.0, 1.0, 0.0])
    e1 = trial - np.dot(trial, normal) * normal
    e1 /= np.linalg.norm(e1)
    e2 = np.cross(normal, e1)
    e2 /= np.linalg.norm(e2)
    return e1, e2


def adsorbate_atoms(kind, n_pos, normal, nh_length):
    """Return symbols and positions for N/NHx with hydrogens pointing outward."""
    symbols = ["N"]
    positions = [np.asarray(n_pos, dtype=float)]
    if kind == "N":
        return symbols, positions

    e1, e2 = make_tangent_basis(normal)

    if kind == "NH":
        directions = [normal]
    elif kind == "NH2":
        # H-N-H angle near 104 deg, symmetric around the outward normal.
        theta = math.radians(52.0)
        directions = [
            math.cos(theta) * normal + math.sin(theta) * e1,
            math.cos(theta) * normal - math.sin(theta) * e1,
        ]
    elif kind == "NH3":
        # Trigonal-pyramidal NH3 with H atoms above the N atom.
        theta = math.radians(68.0)
        directions = []
        for phi in (0.0, 2.0 * math.pi / 3.0, 4.0 * math.pi / 3.0):
            tangent = math.cos(phi) * e1 + math.sin(phi) * e2
            directions.append(math.cos(theta) * normal + math.sin(theta) * tangent)
    else:
        sys.exit(f"[error] Unsupported adsorbate: {kind}")

    for direction in directions:
        direction = np.asarray(direction, dtype=float)
        direction /= np.linalg.norm(direction)
        symbols.append("H")
        positions.append(positions[0] + nh_length * direction)
    return symbols, positions


def add_adsorbate(atoms, site, kind, normal, height, nh_length):
    from ase import Atoms

    out = atoms.copy()
    n_pos = site["plane_pos"] + normal * float(height)
    symbols, positions = adsorbate_atoms(kind, n_pos, normal, nh_length)
    ads = Atoms(symbols=symbols, positions=positions, cell=out.cell, pbc=out.pbc)
    out += ads
    return out


def output_name(input_path, adsorbate, site_index):
    stem, _ = os.path.splitext(input_path)
    return f"{stem}_{adsorbate}_fcc{site_index}.cif"


def print_sites(sites, layers, normal):
    print(f"Detected slab layers: {len(layers)}")
    print(f"Top layer atoms used for surface plane: {len(layers[-1]['indices'])}")
    print(f"Third layer atoms used for FCC hollows: {len(layers[-3]['indices'])}")
    print("FCC sites:")
    for i, site in enumerate(sites, start=1):
        h = float(site["plane_pos"] @ normal)
        fx, fy = site["frac_xy"]
        ax, ay = site["anchor_frac_xy"]
        top_atoms = ",".join(str(idx + 1) for idx in site["top_atoms"]) or "-"
        print(
            f"  {i:3d}: frac_xy=({fx:.10f}, {fy:.10f})  "
            f"anchor=({ax:.10f}, {ay:.10f})  plane_h={h:.6f} A  "
            f"top_atoms={top_atoms}"
        )


def parse_site_selection(text, n_sites):
    if text.lower() == "all":
        return list(range(1, n_sites + 1))
    try:
        idx = int(text)
    except ValueError:
        sys.exit("[error] --site must be an integer index or 'all'.")
    if idx < 1 or idx > n_sites:
        sys.exit(f"[error] --site {idx} out of range; valid range is 1..{n_sites}.")
    return [idx]


def main(argv=None):
    args = parse_args(argv)
    atoms = read_cif(args.cif)
    atoms.set_pbc(True)

    normal = normal_from_cell(atoms.cell, args.side)
    metal_symbols = set(args.metals)
    sites, layers = detect_fcc_sites(
        atoms, normal, args.side, metal_symbols, args.layer_tol, args.site_mode
    )
    if not sites:
        sys.exit("[error] No FCC sites detected.")

    print_sites(sites, layers, normal)
    if args.list_sites:
        return 0

    selected = parse_site_selection(args.site, len(sites))
    if args.output and len(selected) > 1:
        print("[warn] --output is ignored with --site all; writing numbered files.")

    for site_number in selected:
        site = sites[site_number - 1]
        out_atoms = add_adsorbate(
            atoms, site, args.adsorbate, normal, args.height, args.nh
        )
        out_path = args.output if args.output and len(selected) == 1 else (
            output_name(args.cif, args.adsorbate, site_number)
        )
        write_cif(out_path, out_atoms)
        print(
            f"[done] Wrote {out_path}  "
            f"(adsorbate={args.adsorbate}, fcc_site={site_number})"
        )

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
