#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Extract surface-atom d-orbital Hamiltonian blocks from an OpenMX .scfout file.

The script uses:
  1. OpenMX .dat input to reconstruct atom order and orbital labels.
  2. metadata/top_atoms.jsonl from hea_dataset.py for the surface atom order.
  3. OpenMX .scfout binary Hks blocks for Hamiltonian matrix elements.

The output is a compressed NumPy archive:

  H_d[spin, surface_i, surface_j, d_i, d_j]
  d_basis_indices[surface_i, d_i]
  surface_atom_ids[surface_i]
  openmx_atom_indices[surface_i]

Only same-cell blocks are exported by default, which is the right convention for
the finite slab cell used in this workflow.
"""

import argparse
import json
import re
import sqlite3
import struct
import sys
from dataclasses import dataclass
from datetime import datetime, timezone
from pathlib import Path

import numpy as np


ANGULAR_LABELS = {
    "s": ["s"],
    "p": ["px", "py", "pz"],
    "d": ["dxy", "dyz", "dz2", "dxz", "dx2-y2"],
    "f": ["f1", "f2", "f3", "f4", "f5", "f6", "f7"],
}


def utc_now():
    return datetime.now(timezone.utc).isoformat(timespec="seconds")


def read_jsonl(path):
    rows = []
    with open(path, "r", encoding="utf-8") as f:
        for line in f:
            line = line.strip()
            if line:
                rows.append(json.loads(line))
    return rows


def strip_comments(line):
    return line.split("#", 1)[0].strip()


def parse_dat(path):
    species_basis = {}
    atoms = []
    in_species = False
    in_atoms = False
    with open(path, "r", encoding="utf-8") as f:
        for raw in f:
            line = strip_comments(raw)
            if not line:
                continue
            if line.startswith("<Definition.of.Atomic.Species"):
                in_species = True
                continue
            if line.startswith("Definition.of.Atomic.Species>"):
                in_species = False
                continue
            if line.startswith("<Atoms.SpeciesAndCoordinates"):
                in_atoms = True
                continue
            if line.startswith("Atoms.SpeciesAndCoordinates>"):
                in_atoms = False
                continue
            if in_species:
                parts = line.split()
                if len(parts) >= 2:
                    species_basis[parts[0]] = parts[1]
                continue
            if in_atoms:
                parts = line.split()
                if len(parts) >= 2:
                    atoms.append({"openmx_atom_index": int(parts[0]), "element": parts[1]})

    if not species_basis:
        raise SystemExit(f"[error] No Definition.of.Atomic.Species block in {path}")
    if not atoms:
        raise SystemExit(f"[error] No Atoms.SpeciesAndCoordinates block in {path}")
    return species_basis, atoms


def orbital_labels_from_basis(basis_name):
    if "-" not in basis_name:
        raise SystemExit(f"[error] Cannot parse OpenMX basis name: {basis_name}")
    contraction = basis_name.split("-", 1)[1]
    tokens = re.findall(r"([spdf])(\d+)", contraction)
    if not tokens:
        raise SystemExit(f"[error] Cannot parse basis contraction: {basis_name}")

    labels = []
    for angular, count_text in tokens:
        count = int(count_text)
        base = ANGULAR_LABELS[angular]
        for zeta in range(1, count + 1):
            for label in base:
                labels.append(f"{label}_z{zeta}")
    return labels


def build_basis(atoms, species_basis):
    basis = []
    offsets = {}
    cursor = 0
    for atom in atoms:
        element = atom["element"]
        if element not in species_basis:
            raise SystemExit(f"[error] Missing species basis for {element}")
        labels = orbital_labels_from_basis(species_basis[element])
        atom_index = atom["openmx_atom_index"]
        offsets[atom_index] = (cursor, cursor + len(labels))
        for local_i, label in enumerate(labels):
            basis.append({
                "basis_index": cursor + local_i,
                "openmx_atom_index": atom_index,
                "element": element,
                "orbital_label": label,
                "orbital_type": label[0],
            })
        cursor += len(labels)
    return basis, offsets


@dataclass
class ScfoutData:
    atomnum: int
    spin_switch: int
    total_orbitals_by_atom: list
    fnan: list
    natn: list
    ncn: list
    hamiltonian: np.ndarray


class BinaryReader:
    def __init__(self, path, endian="<"):
        self.path = Path(path)
        self.data = self.path.read_bytes()
        self.offset = 0
        self.endian = endian

    def remaining(self):
        return len(self.data) - self.offset

    def read_int(self):
        if self.remaining() < 4:
            raise EOFError("Unexpected end of file while reading int")
        value = struct.unpack_from(self.endian + "i", self.data, self.offset)[0]
        self.offset += 4
        return value

    def read_double(self):
        if self.remaining() < 8:
            raise EOFError("Unexpected end of file while reading double")
        value = struct.unpack_from(self.endian + "d", self.data, self.offset)[0]
        self.offset += 8
        return value

    def read_ints(self, n):
        if n <= 0:
            return []
        nbytes = 4 * n
        if self.remaining() < nbytes:
            raise EOFError("Unexpected end of file while reading int array")
        values = struct.unpack_from(self.endian + f"{n}i", self.data, self.offset)
        self.offset += nbytes
        return list(values)

    def read_doubles(self, n):
        if n <= 0:
            return []
        nbytes = 8 * n
        if self.remaining() < nbytes:
            raise EOFError("Unexpected end of file while reading double array")
        values = struct.unpack_from(self.endian + f"{n}d", self.data, self.offset)
        self.offset += nbytes
        return list(values)


def spin_channels(spin_switch):
    if spin_switch == 0:
        return 1
    if spin_switch == 1:
        return 2
    if spin_switch == 3:
        return 4
    raise ValueError(f"Unsupported OpenMX SpinP_switch={spin_switch}")


def parse_scfout_binary(path, expected_total_orbitals=None):
    errors = []
    for endian in ("<", ">"):
        for vector_width in (4,):
            try:
                return _parse_scfout_hsfileout(
                    path, endian, vector_width, expected_total_orbitals
                )
            except Exception as exc:
                errors.append(
                    f"hsfileout endian={endian} vector_width={vector_width}: {exc}"
                )
        for vector_width in (4, 3):
            try:
                return _parse_scfout_legacy_guess(
                    path, endian, vector_width, expected_total_orbitals
                )
            except Exception as exc:
                errors.append(f"legacy endian={endian} vector_width={vector_width}: {exc}")
    joined = "\n  ".join(errors)
    raise SystemExit(
        "[error] Could not parse .scfout as OpenMX binary layout.\n"
        "Tried common OpenMX 3.9/4.x layout variants:\n  "
        f"{joined}"
    )


def _read_sparse_matrix(r, atomnum, total_orbitals_by_atom, fnan, natn, ncn):
    offsets = np.cumsum([0] + total_orbitals_by_atom)
    total = int(offsets[-1])
    mat = np.zeros((total, total), dtype=np.float64)

    for ct in range(atomnum):
        row0, row1 = int(offsets[ct]), int(offsets[ct + 1])
        tno1 = row1 - row0
        for neigh_i in range(fnan[ct] + 1):
            gh_an = natn[ct][neigh_i]
            cell_i = ncn[ct][neigh_i]
            if gh_an < 1 or gh_an > atomnum:
                raise ValueError(f"neighbor atom out of range: {gh_an}")
            col0, col1 = int(offsets[gh_an - 1]), int(offsets[gh_an])
            tno2 = col1 - col0
            block = np.array(r.read_doubles(tno1 * tno2), dtype=np.float64)
            block = block.reshape((tno1, tno2))
            if cell_i == 0:
                mat[row0:row1, col0:col1] = block
    return mat


def _parse_scfout_hsfileout(path, endian, vector_width, expected_total_orbitals):
    """Parse the OpenMX HS.fileout binary layout.

    OpenMX writes several 1-based C arrays before the sparse matrix blocks:

        atomnum, Matomnum, Catomnum, Latomnum, Ratomnum, TCpyCell, SpinP_switch
        atv, atv_ijk
        Total_NumOrbs, FNAN, natn, ncn
        tv, rtv, Gxyz
        Hks[spin], OLP, ...

    The Hamiltonian blocks are the first matrix blocks after coordinates; for
    spin-polarized calculations the first two blocks are up/down Hks.
    """
    r = BinaryReader(path, endian=endian)
    atomnum = r.read_int()
    matomnum = r.read_int()
    catomnum = r.read_int()
    latomnum = r.read_int()
    ratomnum = r.read_int()
    tcpy_cell = r.read_int()
    spin_switch = r.read_int()

    if atomnum <= 0 or atomnum > 100000:
        raise ValueError(f"implausible atomnum={atomnum}")
    if matomnum <= 0 or matomnum > atomnum:
        raise ValueError(f"implausible Matomnum={matomnum}")
    if min(catomnum, latomnum, ratomnum) < 0:
        raise ValueError("negative atom counts in header")
    if tcpy_cell < 0 or tcpy_cell > 1000000:
        raise ValueError(f"implausible TCpyCell={tcpy_cell}")
    if spin_switch not in (0, 1, 3):
        raise ValueError(f"unsupported or implausible SpinP_switch={spin_switch}")

    r.read_doubles((tcpy_cell + 1) * vector_width)
    r.read_ints((tcpy_cell + 1) * vector_width)

    total_orbitals_by_atom = r.read_ints(atomnum)
    if any(v <= 0 or v > 1000 for v in total_orbitals_by_atom):
        raise ValueError("implausible orbital counts")
    if expected_total_orbitals is not None:
        got = sum(total_orbitals_by_atom)
        if got != expected_total_orbitals:
            raise ValueError(
                f"orbital count mismatch: scfout={got}, dat={expected_total_orbitals}"
            )

    fnan = r.read_ints(atomnum)
    if any(v < 0 or v > atomnum * max(1, tcpy_cell + 1) for v in fnan):
        raise ValueError("implausible FNAN values")

    natn = []
    ncn = []
    for ct in range(atomnum):
        natn.append(r.read_ints(fnan[ct] + 1))
    for ct in range(atomnum):
        ncn.append(r.read_ints(fnan[ct] + 1))

    # tv[1..3][0..3], rtv[1..3][0..3], Gxyz[1..atomnum][0..3]
    r.read_doubles(3 * vector_width)
    r.read_doubles(3 * vector_width)
    r.read_doubles(atomnum * vector_width)

    nspin = spin_channels(spin_switch)
    h_mats = []
    for _ in range(nspin):
        h_mats.append(
            _read_sparse_matrix(
                r, atomnum, total_orbitals_by_atom, fnan, natn, ncn
            )
        )
    h = np.stack(h_mats, axis=0)

    return ScfoutData(
        atomnum=atomnum,
        spin_switch=spin_switch,
        total_orbitals_by_atom=total_orbitals_by_atom,
        fnan=fnan,
        natn=natn,
        ncn=ncn,
        hamiltonian=h,
    )


def _parse_scfout_legacy_guess(path, endian, vector_width, expected_total_orbitals):
    r = BinaryReader(path, endian=endian)
    atomnum = r.read_int()
    spin_switch = r.read_int()
    catomnum = r.read_int()
    latomnum = r.read_int()
    ratomnum = r.read_int()
    tcpy_cell = r.read_int()
    _order_max = r.read_int()

    if atomnum <= 0 or atomnum > 100000:
        raise ValueError(f"implausible atomnum={atomnum}")
    if spin_switch not in (0, 1, 3):
        raise ValueError(f"unsupported or implausible SpinP_switch={spin_switch}")
    if tcpy_cell < 0 or tcpy_cell > 1000000:
        raise ValueError(f"implausible TCpyCell={tcpy_cell}")
    if min(catomnum, latomnum, ratomnum) < 0:
        raise ValueError("negative atom counts in header")

    _what_species = r.read_ints(atomnum)
    total_orbitals_by_atom = r.read_ints(atomnum)
    if any(v <= 0 or v > 1000 for v in total_orbitals_by_atom):
        raise ValueError("implausible orbital counts")
    if expected_total_orbitals is not None:
        got = sum(total_orbitals_by_atom)
        if got != expected_total_orbitals:
            raise ValueError(
                f"orbital count mismatch: scfout={got}, dat={expected_total_orbitals}"
            )

    fnan = r.read_ints(atomnum)
    if any(v < 0 or v > atomnum * max(1, tcpy_cell + 1) for v in fnan):
        raise ValueError("implausible FNAN values")

    natn = []
    ncn = []
    for ct in range(atomnum):
        natn.append(r.read_ints(fnan[ct] + 1))
    for ct in range(atomnum):
        ncn.append(r.read_ints(fnan[ct] + 1))

    # Translation vectors and cell matrices are present before Hks. Some OpenMX
    # builds write four components because the original C arrays are 1-based.
    r.read_doubles((tcpy_cell + 1) * vector_width)
    r.read_ints((tcpy_cell + 1) * vector_width)
    r.read_doubles(3 * vector_width)
    r.read_doubles(3 * vector_width)
    r.read_doubles(atomnum * vector_width)

    nspin = spin_channels(spin_switch)
    h_mats = []
    for _ in range(nspin):
        h_mats.append(
            _read_sparse_matrix(
                r, atomnum, total_orbitals_by_atom, fnan, natn, ncn
            )
        )
    h = np.stack(h_mats, axis=0)

    return ScfoutData(
        atomnum=atomnum,
        spin_switch=spin_switch,
        total_orbitals_by_atom=total_orbitals_by_atom,
        fnan=fnan,
        natn=natn,
        ncn=ncn,
        hamiltonian=h,
    )


def d_indices_for_surface(top_rows, basis, offsets):
    basis_by_index = {row["basis_index"]: row for row in basis}
    atom_ids = []
    openmx_indices = []
    d_lists = []
    d_label_lists = []

    ordered = sorted(top_rows, key=lambda row: (row["row"], row["col"]))
    for row in ordered:
        atom_index = row.get("openmx_atom_index")
        if atom_index is None:
            atom_index = row.get("relaxed_ase_index", row["initial_ase_index"]) + 1
        atom_index = int(atom_index)
        if atom_index not in offsets:
            raise SystemExit(f"[error] OpenMX atom index not found: {atom_index}")

        start, end = offsets[atom_index]
        indices = []
        labels = []
        for basis_index in range(start, end):
            brow = basis_by_index[basis_index]
            if brow["orbital_type"] == "d":
                indices.append(basis_index)
                labels.append(brow["orbital_label"])
        atom_ids.append(row["atom_id"])
        openmx_indices.append(atom_index)
        d_lists.append(indices)
        d_label_lists.append(labels)

    return atom_ids, openmx_indices, d_lists, d_label_lists


def extract_blocks(h, d_lists):
    nspin = h.shape[0]
    nsurf = len(d_lists)
    max_d = max((len(items) for items in d_lists), default=0)
    if max_d == 0:
        raise SystemExit("[error] No d orbitals found for the selected surface atoms.")

    out = np.full((nspin, nsurf, nsurf, max_d, max_d), np.nan, dtype=np.float64)
    d_basis = np.full((nsurf, max_d), -1, dtype=np.int64)
    for i, rows in enumerate(d_lists):
        d_basis[i, :len(rows)] = rows
        for j, cols in enumerate(d_lists):
            if rows and cols:
                out[:, i, j, :len(rows), :len(cols)] = h[:, rows][:, :, cols]
    return out, d_basis


def write_basis_jsonl(path, basis):
    path.parent.mkdir(parents=True, exist_ok=True)
    with open(path, "w", encoding="utf-8") as f:
        for row in basis:
            f.write(json.dumps(row, sort_keys=True, separators=(",", ":")) + "\n")


def register_export(args, n_surface_atoms, spin_channels_count):
    db = Path(args.dataset_root) / "index.sqlite"
    if not db.is_file() or not args.surface_id:
        return
    with sqlite3.connect(db) as con:
        con.execute(
            """
            INSERT INTO hamiltonian_exports(
                surface_id, scfout_path, output_npz, n_surface_atoms,
                spin_channels, created_at
            ) VALUES (?, ?, ?, ?, ?, ?)
            """,
            (
                args.surface_id, str(args.scfout), str(args.output),
                n_surface_atoms, spin_channels_count, utc_now(),
            ),
        )
        con.commit()


def main(argv=None):
    p = argparse.ArgumentParser(
        description="Extract surface d-orbital Hamiltonian blocks from OpenMX .scfout.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    p.add_argument("--scfout", required=True, help="OpenMX .scfout file.")
    p.add_argument("--dat", required=True, help="OpenMX .dat input used for the run.")
    p.add_argument(
        "--top-atoms", required=True,
        help="metadata/top_atoms.jsonl generated by hea_dataset.py.",
    )
    p.add_argument(
        "-o", "--output", required=True,
        help="Output compressed NumPy archive, e.g. hamiltonian_d_surface.npz.",
    )
    p.add_argument(
        "--basis-output", default=None,
        help="Optional JSONL basis mapping output. Defaults to <output>.basis.jsonl.",
    )
    p.add_argument("--dataset-root", default="dataset")
    p.add_argument("--surface-id", default=None)
    args = p.parse_args(argv)

    species_basis, atoms = parse_dat(args.dat)
    basis, offsets = build_basis(atoms, species_basis)
    expected_total = len(basis)
    scfout = parse_scfout_binary(args.scfout, expected_total_orbitals=expected_total)

    top_rows = read_jsonl(args.top_atoms)
    atom_ids, openmx_indices, d_lists, d_label_lists = d_indices_for_surface(
        top_rows, basis, offsets
    )
    h_d, d_basis = extract_blocks(scfout.hamiltonian, d_lists)

    output = Path(args.output)
    output.parent.mkdir(parents=True, exist_ok=True)
    max_d = d_basis.shape[1]
    d_labels = np.full((len(d_label_lists), max_d), "", dtype="<U32")
    for i, labels in enumerate(d_label_lists):
        d_labels[i, :len(labels)] = labels

    np.savez_compressed(
        output,
        H_d=h_d,
        d_basis_indices=d_basis,
        d_labels=d_labels,
        surface_atom_ids=np.array(atom_ids, dtype="<U64"),
        openmx_atom_indices=np.array(openmx_indices, dtype=np.int64),
        spin_switch=np.array([scfout.spin_switch], dtype=np.int64),
        source_scfout=np.array([str(args.scfout)], dtype="<U1024"),
        source_dat=np.array([str(args.dat)], dtype="<U1024"),
    )

    basis_output = (
        Path(args.basis_output)
        if args.basis_output
        else output.with_suffix(output.suffix + ".basis.jsonl")
    )
    write_basis_jsonl(basis_output, basis)
    register_export(args, len(atom_ids), h_d.shape[0])

    print(f"[done] Wrote {output}")
    print(f"[done] Wrote {basis_output}")
    print(
        f"[info] surface_atoms={len(atom_ids)} spin_channels={h_d.shape[0]} "
        f"d_orbitals_per_atom(max)={max_d}"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
