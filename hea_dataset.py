#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Dataset manager for HEA surface training data.

The dataset identity follows the convention requested for this project:

    surface_id = element names + element fractions only
    sample_id  = one concrete structure under that composition

All structural files are kept on disk, while the searchable state lives in a
SQLite index. Per-sample metadata is also written as JSONL/NPY files so that
training code can load it without going through SQL.
"""

import argparse
import hashlib
import json
import os
import shutil
import sqlite3
import sys
from datetime import datetime, timezone
from pathlib import Path

import numpy as np


SCHEMA_VERSION = 1
METALS_DEFAULT = [
    "Fe", "Co", "Ni", "Cu", "Mo", "Zn", "Ga", "In", "Sn", "W",
    "Cr", "Mn", "Pd", "Pt", "Rh", "Ir",
]


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


def write_json(path, data):
    path.parent.mkdir(parents=True, exist_ok=True)
    with open(path, "w", encoding="utf-8") as f:
        json.dump(data, f, indent=2, sort_keys=True)
        f.write("\n")


def write_jsonl(path, rows):
    path.parent.mkdir(parents=True, exist_ok=True)
    with open(path, "w", encoding="utf-8") as f:
        for row in rows:
            f.write(json.dumps(row, sort_keys=True, separators=(",", ":")) + "\n")


def sha256_file(path):
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for chunk in iter(lambda: f.read(1024 * 1024), b""):
            h.update(chunk)
    return h.hexdigest()


def parse_composition(items):
    pairs = []
    for item in items:
        if "=" not in item:
            raise SystemExit(
                f"[error] Bad composition entry '{item}', expected EL=FRACTION."
            )
        element, value = item.split("=", 1)
        element = element.strip()
        if not element:
            raise SystemExit(f"[error] Bad composition entry '{item}'.")
        try:
            fraction = float(value)
        except ValueError:
            raise SystemExit(f"[error] Bad fraction in '{item}'.")
        if fraction < 0:
            raise SystemExit(f"[error] Fraction must be non-negative in '{item}'.")
        pairs.append((element, fraction))

    total = sum(v for _, v in pairs)
    if total <= 0:
        raise SystemExit("[error] Composition fractions must sum to > 0.")
    return [(el, 100.0 * value / total) for el, value in pairs]


def format_fraction(value):
    text = f"{value:.4f}".rstrip("0").rstrip(".")
    return text if text else "0"


def surface_id_from_composition(composition):
    # Preserve user-specified element order; it usually follows the alloy name.
    return "-".join(f"{el}_{format_fraction(frac)}" for el, frac in composition)


def dataset_paths(root, surface_id=None, sample_id=None):
    root = Path(root)
    paths = {
        "root": root,
        "db": root / "index.sqlite",
        "manifest": root / "dataset_manifest.json",
    }
    if surface_id:
        paths["surface"] = root / "surfaces" / surface_id
    if surface_id and sample_id:
        sample = root / "surfaces" / surface_id / sample_id
        paths.update({
            "sample": sample,
            "structures": sample / "structures",
            "metadata": sample / "metadata",
            "openmx_slab": sample / "openmx_slab",
            "adsorbates": sample / "adsorbates",
            "sample_manifest": sample / "manifest.json",
        })
    return paths


def connect_db(root):
    return sqlite3.connect(dataset_paths(root)["db"])


def init_db(root):
    paths = dataset_paths(root)
    paths["root"].mkdir(parents=True, exist_ok=True)
    with sqlite3.connect(paths["db"]) as con:
        con.executescript(
            """
            PRAGMA foreign_keys = ON;

            CREATE TABLE IF NOT EXISTS dataset_info (
                key TEXT PRIMARY KEY,
                value TEXT NOT NULL
            );

            CREATE TABLE IF NOT EXISTS surfaces (
                surface_id TEXT PRIMARY KEY,
                composition_json TEXT NOT NULL,
                created_at TEXT NOT NULL
            );

            CREATE TABLE IF NOT EXISTS samples (
                surface_id TEXT NOT NULL,
                sample_id TEXT NOT NULL,
                path TEXT NOT NULL,
                status TEXT NOT NULL DEFAULT 'created',
                initial_cif TEXT,
                relaxed_cif TEXT,
                created_at TEXT NOT NULL,
                updated_at TEXT NOT NULL,
                PRIMARY KEY (surface_id, sample_id),
                FOREIGN KEY (surface_id) REFERENCES surfaces(surface_id)
            );

            CREATE TABLE IF NOT EXISTS artifacts (
                artifact_id INTEGER PRIMARY KEY AUTOINCREMENT,
                surface_id TEXT NOT NULL,
                sample_id TEXT NOT NULL,
                kind TEXT NOT NULL,
                path TEXT NOT NULL,
                sha256 TEXT,
                created_at TEXT NOT NULL,
                FOREIGN KEY (surface_id, sample_id)
                    REFERENCES samples(surface_id, sample_id)
            );

            CREATE TABLE IF NOT EXISTS top_atoms (
                surface_id TEXT NOT NULL,
                sample_id TEXT NOT NULL,
                atom_id TEXT NOT NULL,
                row INTEGER NOT NULL,
                col INTEGER NOT NULL,
                element TEXT NOT NULL,
                initial_ase_index INTEGER NOT NULL,
                relaxed_ase_index INTEGER,
                openmx_atom_index INTEGER,
                initial_x REAL,
                initial_y REAL,
                initial_z REAL,
                relaxed_x REAL,
                relaxed_y REAL,
                relaxed_z REAL,
                PRIMARY KEY (surface_id, sample_id, atom_id),
                FOREIGN KEY (surface_id, sample_id)
                    REFERENCES samples(surface_id, sample_id)
            );

            CREATE TABLE IF NOT EXISTS fcc_sites (
                surface_id TEXT NOT NULL,
                sample_id TEXT NOT NULL,
                site_id TEXT NOT NULL,
                site_index INTEGER NOT NULL,
                row INTEGER NOT NULL,
                col INTEGER NOT NULL,
                site_type TEXT NOT NULL DEFAULT 'fcc',
                frac_x REAL NOT NULL,
                frac_y REAL NOT NULL,
                plane_x REAL NOT NULL,
                plane_y REAL NOT NULL,
                plane_z REAL NOT NULL,
                top_atom_ids_json TEXT NOT NULL,
                adsorption_energy_eV REAL,
                energy_status TEXT NOT NULL DEFAULT 'empty',
                PRIMARY KEY (surface_id, sample_id, site_id),
                FOREIGN KEY (surface_id, sample_id)
                    REFERENCES samples(surface_id, sample_id)
            );

            CREATE TABLE IF NOT EXISTS adsorbate_configs (
                config_id TEXT PRIMARY KEY,
                surface_id TEXT NOT NULL,
                sample_id TEXT NOT NULL,
                site_id TEXT NOT NULL,
                adsorbate TEXT NOT NULL,
                path TEXT NOT NULL,
                initial_cif TEXT,
                relaxed_cif TEXT,
                adsorption_energy_eV REAL,
                energy_status TEXT NOT NULL DEFAULT 'empty',
                updated_at TEXT NOT NULL,
                FOREIGN KEY (surface_id, sample_id, site_id)
                    REFERENCES fcc_sites(surface_id, sample_id, site_id)
            );

            CREATE TABLE IF NOT EXISTS hamiltonian_exports (
                export_id INTEGER PRIMARY KEY AUTOINCREMENT,
                surface_id TEXT NOT NULL,
                sample_id TEXT NOT NULL,
                scfout_path TEXT NOT NULL,
                output_npz TEXT NOT NULL,
                n_surface_atoms INTEGER NOT NULL,
                spin_channels INTEGER NOT NULL,
                created_at TEXT NOT NULL,
                FOREIGN KEY (surface_id, sample_id)
                    REFERENCES samples(surface_id, sample_id)
            );
            """
        )
        con.execute(
            "INSERT OR REPLACE INTO dataset_info(key, value) VALUES (?, ?)",
            ("schema_version", str(SCHEMA_VERSION)),
        )
        con.commit()

    write_json(paths["manifest"], {
        "schema_version": SCHEMA_VERSION,
        "created_or_updated_at": utc_now(),
        "layout": "surfaces/<surface_id>/<sample_id>",
        "surface_id_rule": "element_fraction_pairs_only",
    })
    print(f"[done] Initialized dataset at {paths['root']}")


def copy_if_requested(src, dst):
    if src is None:
        return None
    src = Path(src)
    if not src.is_file():
        raise SystemExit(f"[error] File not found: {src}")
    dst.parent.mkdir(parents=True, exist_ok=True)
    shutil.copy2(src, dst)
    return dst


def create_sample(args):
    init_db(args.root)
    composition = parse_composition(args.composition)
    surface_id = args.surface_id or surface_id_from_composition(composition)
    sample_id = args.sample_id

    paths = dataset_paths(args.root, surface_id, sample_id)
    for key in ("structures", "metadata", "openmx_slab", "adsorbates"):
        paths[key].mkdir(parents=True, exist_ok=True)

    initial = copy_if_requested(
        args.initial_cif, paths["structures"] / "00_initial_sqs.cif"
    )
    relaxed = copy_if_requested(
        args.relaxed_cif, paths["structures"] / "01_relaxed_slab.cif"
    )

    manifest = {
        "surface_id": surface_id,
        "sample_id": sample_id,
        "composition": [{"element": el, "fraction_percent": frac}
                        for el, frac in composition],
        "created_at": utc_now(),
        "files": {
            "initial_cif": relpath_or_none(initial, paths["sample"]),
            "relaxed_cif": relpath_or_none(relaxed, paths["sample"]),
            "top_atoms": "metadata/top_atoms.jsonl",
            "fcc_sites": "metadata/fcc_sites.jsonl",
            "atom_grid": "metadata/atom_grid.npy",
            "site_grid": "metadata/site_grid.npy",
        },
    }
    write_json(paths["sample_manifest"], manifest)

    with connect_db(args.root) as con:
        con.execute(
            "INSERT OR IGNORE INTO surfaces(surface_id, composition_json, created_at) "
            "VALUES (?, ?, ?)",
            (surface_id, json.dumps(manifest["composition"]), utc_now()),
        )
        con.execute(
            """
            INSERT OR REPLACE INTO samples(
                surface_id, sample_id, path, status, initial_cif, relaxed_cif,
                created_at, updated_at
            ) VALUES (?, ?, ?, ?, ?, ?, ?, ?)
            """,
            (
                surface_id, sample_id, str(paths["sample"]),
                "created", relpath_or_none(initial, paths["sample"]),
                relpath_or_none(relaxed, paths["sample"]), utc_now(), utc_now(),
            ),
        )
        for kind, artifact in (("initial_cif", initial), ("relaxed_cif", relaxed)):
            if artifact is not None:
                con.execute(
                    """
                    INSERT INTO artifacts(surface_id, sample_id, kind, path, sha256, created_at)
                    VALUES (?, ?, ?, ?, ?, ?)
                    """,
                    (
                        surface_id, sample_id, kind,
                        str(artifact), sha256_file(artifact), utc_now(),
                    ),
                )
        con.commit()

    print(f"[done] Created sample {surface_id}/{sample_id}")


def relpath_or_none(path, base):
    if path is None:
        return None
    return os.path.relpath(path, base).replace("\\", "/")


def read_cif(path):
    try:
        from ase.io import read
    except ImportError:
        raise SystemExit("[error] ASE is required for CIF indexing: pip install ase")
    if not Path(path).is_file():
        raise SystemExit(f"[error] CIF file not found: {path}")
    atoms = read(path)
    if atoms.cell.rank != 3:
        raise SystemExit("[error] The CIF must contain a full 3D cell.")
    return atoms


def normal_from_cell(cell, side):
    normal = np.asarray(cell[2], dtype=float)
    norm = np.linalg.norm(normal)
    if norm < 1.0e-12:
        raise SystemExit("[error] Cell c vector is zero; cannot define surface normal.")
    normal = normal / norm
    return -normal if side == "bottom" else normal


def group_layers(atoms, normal, metal_symbols, tol):
    symbols = np.array(atoms.get_chemical_symbols())
    positions = atoms.get_positions()
    metal_mask = np.array([sym in metal_symbols for sym in symbols])
    if not np.any(metal_mask):
        raise SystemExit("[error] No slab atoms found. Check --metals.")

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
    return layers


def assign_grid(indices, atoms):
    scaled = atoms.get_scaled_positions(wrap=True)
    ordered = sorted(
        indices,
        key=lambda idx: (
            round(float(scaled[idx, 1]), 10),
            round(float(scaled[idx, 0]), 10),
        ),
    )

    rows = []
    for idx in ordered:
        fy = float(scaled[idx, 1])
        if not rows or abs(fy - rows[-1]["fy"]) > 1.0e-5:
            rows.append({"fy": fy, "indices": [idx]})
        else:
            rows[-1]["indices"].append(idx)

    mapping = {}
    max_cols = 0
    for row_idx, row in enumerate(rows):
        row["indices"].sort(key=lambda idx: float(scaled[idx, 0]))
        max_cols = max(max_cols, len(row["indices"]))
        for col_idx, atom_idx in enumerate(row["indices"]):
            mapping[int(atom_idx)] = (row_idx, col_idx)
    return mapping, len(rows), max_cols


def assign_xy_grid(items, xy_getter):
    ordered = sorted(
        range(len(items)),
        key=lambda i: (
            round(float(xy_getter(items[i])[1]), 10),
            round(float(xy_getter(items[i])[0]), 10),
        ),
    )
    rows = []
    for item_i in ordered:
        _, fy = xy_getter(items[item_i])
        fy = float(fy)
        if not rows or abs(fy - rows[-1]["fy"]) > 1.0e-5:
            rows.append({"fy": fy, "indices": [item_i]})
        else:
            rows[-1]["indices"].append(item_i)

    mapping = {}
    max_cols = 0
    for row_idx, row in enumerate(rows):
        row["indices"].sort(key=lambda i: float(xy_getter(items[i])[0]))
        max_cols = max(max_cols, len(row["indices"]))
        for col_idx, item_i in enumerate(row["indices"]):
            mapping[item_i] = (row_idx, col_idx)
    return mapping, len(rows), max_cols


def index_surface(args):
    paths = dataset_paths(args.root, args.surface_id, args.sample_id)
    initial_cif = Path(args.initial_cif or paths["structures"] / "00_initial_sqs.cif")
    relaxed_cif = Path(args.relaxed_cif or paths["structures"] / "01_relaxed_slab.cif")
    initial_atoms = read_cif(initial_cif)
    relaxed_atoms = read_cif(relaxed_cif) if relaxed_cif.is_file() else None

    normal = normal_from_cell(initial_atoms.cell, args.side)
    layers = group_layers(initial_atoms, normal, set(args.metals), args.layer_tol)
    if not layers:
        raise SystemExit("[error] No layers detected.")

    top_indices = layers[-1]["indices"]
    grid, n_rows, n_cols = assign_grid(top_indices, initial_atoms)
    symbols = initial_atoms.get_chemical_symbols()
    initial_pos = initial_atoms.get_positions()
    relaxed_pos = relaxed_atoms.get_positions() if relaxed_atoms is not None else None
    initial_scaled = initial_atoms.get_scaled_positions(wrap=True)
    relaxed_scaled = (
        relaxed_atoms.get_scaled_positions(wrap=True)
        if relaxed_atoms is not None else None
    )

    rows = []
    atom_grid = np.empty((n_rows, n_cols), dtype="<U32")
    atom_grid[:, :] = ""
    for atom_idx in sorted(top_indices, key=lambda idx: grid[idx]):
        row, col = grid[atom_idx]
        atom_id = f"A_{row:03d}_{col:03d}"
        atom_grid[row, col] = atom_id
        item = {
            "surface_id": args.surface_id,
            "sample_id": args.sample_id,
            "atom_id": atom_id,
            "row": row,
            "col": col,
            "element": symbols[atom_idx],
            "layer_index": len(layers) - 1,
            "initial_ase_index": int(atom_idx),
            "relaxed_ase_index": int(atom_idx) if relaxed_atoms is not None else None,
            "openmx_atom_index": int(atom_idx) + 1,
            "initial_frac_x": float(initial_scaled[atom_idx, 0]),
            "initial_frac_y": float(initial_scaled[atom_idx, 1]),
            "initial_frac_z": float(initial_scaled[atom_idx, 2]),
            "initial_x": float(initial_pos[atom_idx, 0]),
            "initial_y": float(initial_pos[atom_idx, 1]),
            "initial_z": float(initial_pos[atom_idx, 2]),
            "relaxed_frac_x": (
                float(relaxed_scaled[atom_idx, 0]) if relaxed_scaled is not None else None
            ),
            "relaxed_frac_y": (
                float(relaxed_scaled[atom_idx, 1]) if relaxed_scaled is not None else None
            ),
            "relaxed_frac_z": (
                float(relaxed_scaled[atom_idx, 2]) if relaxed_scaled is not None else None
            ),
            "relaxed_x": float(relaxed_pos[atom_idx, 0]) if relaxed_pos is not None else None,
            "relaxed_y": float(relaxed_pos[atom_idx, 1]) if relaxed_pos is not None else None,
            "relaxed_z": float(relaxed_pos[atom_idx, 2]) if relaxed_pos is not None else None,
        }
        rows.append(item)

    top_atoms_path = paths["metadata"] / "top_atoms.jsonl"
    write_jsonl(top_atoms_path, rows)
    np.save(paths["metadata"] / "atom_grid.npy", atom_grid)

    with connect_db(args.root) as con:
        con.execute(
            "DELETE FROM top_atoms WHERE surface_id=? AND sample_id=?",
            (args.surface_id, args.sample_id),
        )
        for item in rows:
            con.execute(
                """
                INSERT INTO top_atoms(
                    surface_id, sample_id, atom_id, row, col, element,
                    initial_ase_index, relaxed_ase_index, openmx_atom_index,
                    initial_x, initial_y, initial_z, relaxed_x, relaxed_y, relaxed_z
                ) VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)
                """,
                (
                    item["surface_id"], item["sample_id"], item["atom_id"],
                    item["row"], item["col"], item["element"],
                    item["initial_ase_index"], item["relaxed_ase_index"],
                    item["openmx_atom_index"], item["initial_x"], item["initial_y"],
                    item["initial_z"], item["relaxed_x"], item["relaxed_y"],
                    item["relaxed_z"],
                ),
            )
        con.execute(
            "UPDATE samples SET status=?, updated_at=? WHERE surface_id=? AND sample_id=?",
            ("surface_indexed", utc_now(), args.surface_id, args.sample_id),
        )
        con.commit()

    print(f"[done] Wrote {top_atoms_path} ({len(rows)} top atoms)")


def detect_sites(args):
    try:
        import add_fcc_adsorbate as ads
    except ImportError as exc:
        raise SystemExit(f"[error] Could not import add_fcc_adsorbate.py: {exc}")

    paths = dataset_paths(args.root, args.surface_id, args.sample_id)
    cif = Path(args.cif or paths["structures"] / "01_relaxed_slab.cif")
    atoms = ads.read_cif(str(cif))
    atoms.set_pbc(True)
    normal = ads.normal_from_cell(atoms.cell, args.side)
    sites, _ = ads.detect_fcc_sites(
        atoms,
        normal,
        args.side,
        set(args.metals),
        args.layer_tol,
        args.site_mode,
    )

    top_atoms_path = paths["metadata"] / "top_atoms.jsonl"
    top_rows = read_jsonl(top_atoms_path) if top_atoms_path.is_file() else []
    by_ase = {row["relaxed_ase_index"]: row for row in top_rows}

    grid, n_rows, n_cols = assign_xy_grid(sites, lambda item: item["frac_xy"])
    rows = []
    site_grid = np.empty((n_rows, n_cols), dtype="<U32")
    site_grid[:, :] = ""
    for i, site in enumerate(sites, start=1):
        row, col = grid[i - 1]
        site_id = f"site_{i:04d}"
        site_grid[row, col] = site_id
        top_atom_ids = [
            by_ase[idx]["atom_id"]
            for idx in site.get("top_atoms", [])
            if idx in by_ase
        ]
        fx, fy = site["frac_xy"]
        px, py, pz = [float(v) for v in site["plane_pos"]]
        rows.append({
            "surface_id": args.surface_id,
            "sample_id": args.sample_id,
            "site_id": site_id,
            "site_index": i,
            "row": row,
            "col": col,
            "site_type": "fcc",
            "anchor_frac_x": float(site["anchor_frac_xy"][0]),
            "anchor_frac_y": float(site["anchor_frac_xy"][1]),
            "frac_x": float(fx),
            "frac_y": float(fy),
            "plane_x": px,
            "plane_y": py,
            "plane_z": pz,
            "source_atom_index": int(site["source_atom"]),
            "top_atom_ids": top_atom_ids,
            "adsorption_energy_eV": None,
            "energy_status": "empty",
        })

    fcc_path = paths["metadata"] / "fcc_sites.jsonl"
    write_jsonl(fcc_path, rows)
    np.save(paths["metadata"] / "site_grid.npy", site_grid)

    with connect_db(args.root) as con:
        con.execute(
            "DELETE FROM fcc_sites WHERE surface_id=? AND sample_id=?",
            (args.surface_id, args.sample_id),
        )
        for item in rows:
            con.execute(
                """
                INSERT INTO fcc_sites(
                    surface_id, sample_id, site_id, site_index, row, col,
                    site_type, frac_x, frac_y, plane_x, plane_y, plane_z,
                    top_atom_ids_json, adsorption_energy_eV, energy_status
                ) VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)
                """,
                (
                    item["surface_id"], item["sample_id"], item["site_id"],
                    item["site_index"], item["row"], item["col"],
                    item["site_type"], item["frac_x"], item["frac_y"],
                    item["plane_x"], item["plane_y"], item["plane_z"],
                    json.dumps(item["top_atom_ids"]),
                    item["adsorption_energy_eV"], item["energy_status"],
                ),
            )
        con.execute(
            "UPDATE samples SET status=?, updated_at=? WHERE surface_id=? AND sample_id=?",
            ("sites_detected", utc_now(), args.surface_id, args.sample_id),
        )
        con.commit()

    print(f"[done] Wrote {fcc_path} ({len(rows)} FCC sites)")


def create_adsorbate_records(args):
    paths = dataset_paths(args.root, args.surface_id, args.sample_id)
    fcc_rows = read_jsonl(paths["metadata"] / "fcc_sites.jsonl")
    ads_root = paths["adsorbates"] / args.adsorbate
    ads_root.mkdir(parents=True, exist_ok=True)

    with connect_db(args.root) as con:
        for site in fcc_rows:
            site_id = site["site_id"]
            config_id = f"{args.surface_id}_{args.sample_id}_{args.adsorbate}_{site_id}"
            config_path = ads_root / site_id
            config_path.mkdir(parents=True, exist_ok=True)
            record = {
                "surface_id": args.surface_id,
                "sample_id": args.sample_id,
                "adsorbate": args.adsorbate,
                "site_id": site_id,
                "config_id": config_id,
                "initial_cif": "00_initial_adsorbate.cif",
                "relaxed_cif": "01_relaxed_adsorbate.cif",
                "adsorption_energy_eV": None,
                "energy_status": "empty",
                "notes": "",
            }
            write_json(config_path / "adsorption_energy.json", record)
            con.execute(
                """
                INSERT OR REPLACE INTO adsorbate_configs(
                    config_id, surface_id, sample_id, site_id, adsorbate, path,
                    initial_cif, relaxed_cif, adsorption_energy_eV, energy_status,
                    updated_at
                ) VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)
                """,
                (
                    config_id, args.surface_id, args.sample_id, site_id,
                    args.adsorbate, str(config_path),
                    "00_initial_adsorbate.cif", "01_relaxed_adsorbate.cif",
                    None, "empty", utc_now(),
                ),
            )
        con.commit()
    print(f"[done] Created adsorbate records under {ads_root}")


def record_energy(args):
    paths = dataset_paths(args.root, args.surface_id, args.sample_id)
    config_path = paths["adsorbates"] / args.adsorbate / args.site_id
    energy_path = config_path / "adsorption_energy.json"
    if not energy_path.is_file():
        raise SystemExit(f"[error] Missing energy record: {energy_path}")
    with open(energy_path, "r", encoding="utf-8") as f:
        record = json.load(f)
    record["adsorption_energy_eV"] = args.energy
    record["energy_status"] = args.status
    record["notes"] = args.notes or record.get("notes", "")
    write_json(energy_path, record)

    with connect_db(args.root) as con:
        con.execute(
            """
            UPDATE adsorbate_configs
            SET adsorption_energy_eV=?, energy_status=?, updated_at=?
            WHERE surface_id=? AND sample_id=? AND adsorbate=? AND site_id=?
            """,
            (
                args.energy, args.status, utc_now(), args.surface_id,
                args.sample_id, args.adsorbate, args.site_id,
            ),
        )
        con.execute(
            """
            UPDATE fcc_sites
            SET adsorption_energy_eV=?, energy_status=?
            WHERE surface_id=? AND sample_id=? AND site_id=?
            """,
            (args.energy, args.status, args.surface_id, args.sample_id, args.site_id),
        )
        con.commit()
    print(f"[done] Recorded {args.energy} eV for {args.adsorbate}/{args.site_id}")


def add_common_sample_args(p):
    p.add_argument("--root", default="dataset", help="Dataset root directory.")
    p.add_argument("--surface-id", required=True, help="Composition-based surface ID.")
    p.add_argument("--sample-id", required=True, help="Concrete sample ID.")


def build_parser():
    p = argparse.ArgumentParser(
        description="Manage HEA surface dataset layout, metadata and SQLite index.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    sub = p.add_subparsers(dest="cmd", required=True)

    sp = sub.add_parser("init", help="Initialize dataset directory and SQLite index.")
    sp.add_argument("--root", default="dataset")
    sp.set_defaults(func=lambda args: init_db(args.root))

    sp = sub.add_parser("create-sample", help="Create a composition/sample directory.")
    sp.add_argument("--root", default="dataset")
    sp.add_argument(
        "--composition", nargs="+", required=True,
        help="Composition entries, e.g. Fe=20 Co=20 Ni=20 Cr=20 Mn=20.",
    )
    sp.add_argument("--surface-id", default=None, help="Override generated surface ID.")
    sp.add_argument("--sample-id", default="sample_0001")
    sp.add_argument("--initial-cif", default=None)
    sp.add_argument("--relaxed-cif", default=None)
    sp.set_defaults(func=create_sample)

    sp = sub.add_parser("index-surface", help="Build top atom metadata/grid.")
    add_common_sample_args(sp)
    sp.add_argument("--initial-cif", default=None)
    sp.add_argument("--relaxed-cif", default=None)
    sp.add_argument("--side", choices=["top", "bottom"], default="top")
    sp.add_argument("--layer-tol", type=float, default=0.60)
    sp.add_argument(
        "--metals", nargs="+",
        default=METALS_DEFAULT,
    )
    sp.set_defaults(func=index_surface)

    sp = sub.add_parser("detect-sites", help="Build FCC site metadata/grid.")
    add_common_sample_args(sp)
    sp.add_argument("--cif", default=None, help="Relaxed clean slab CIF.")
    sp.add_argument("--side", choices=["top", "bottom"], default="top")
    sp.add_argument("--layer-tol", type=float, default=0.60)
    sp.add_argument("--site-mode", choices=["relaxed", "projected"], default="relaxed")
    sp.add_argument(
        "--metals", nargs="+",
        default=METALS_DEFAULT,
    )
    sp.set_defaults(func=detect_sites)

    sp = sub.add_parser(
        "create-adsorbate-records",
        help="Create per-site adsorbate folders and empty energy JSON files.",
    )
    add_common_sample_args(sp)
    sp.add_argument("--adsorbate", required=True, help="Adsorbate name, e.g. N.")
    sp.set_defaults(func=create_adsorbate_records)

    sp = sub.add_parser("record-energy", help="Record a manual adsorption energy.")
    add_common_sample_args(sp)
    sp.add_argument("--adsorbate", required=True)
    sp.add_argument("--site-id", required=True)
    sp.add_argument("--energy", type=float, required=True)
    sp.add_argument("--status", default="manually_entered")
    sp.add_argument("--notes", default="")
    sp.set_defaults(func=record_energy)

    return p


def main(argv=None):
    args = build_parser().parse_args(argv)
    args.func(args)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
