#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Dataset manager for HEA surface training data.

The dataset identity follows the convention requested for this project:

    surface_id = element names + atom counts only
    surface_id = one concrete structure for that composition

All structural files are kept on disk, while the searchable state lives in a
SQLite index. Per-surface metadata is also written as JSONL/NPY files so that
training code can load it without going through SQL.
"""

import argparse
from collections import Counter
import hashlib
import json
import os
import shutil
import sqlite3
import sys
from datetime import datetime, timezone
from pathlib import Path

import numpy as np


SCHEMA_VERSION = 2
DEFAULT_MAX_ATTEMPTS = 500
METALS_DEFAULT = [
    "Fe", "Co", "Ni", "Cu", "Mo", "Zn", "Ga", "In", "Sn", "W",
    "Cr", "Mn", "Pd", "Pt", "Rh", "Ir",
]

# ---------------------------------------------------------------------------
# HEA thermodynamic screening  (Zhang et al. criteria)
# ---------------------------------------------------------------------------

# Universal gas constant [J mol⁻¹ K⁻¹]
_R = 8.314

# Screening thresholds:
_HEA_ENTROPY_MIN    =  1.5 * _R  # ΔS_mix > 1.5 R  [J mol⁻¹ K⁻¹]
_HEA_SIZE_DELTA_MAX =  6.6       # δ ≤ 6.6 %
_HEA_HMIX_MIN       = -15.0      # ΔH_mix lower bound  [kJ mol⁻¹]
_HEA_HMIX_MAX       =  5.0       # ΔH_mix upper bound  [kJ mol⁻¹]
_HEA_VEC_MIN        =  8.0       # VEC ≥ 8  (fcc-stabilizing regime)

# Metallic atomic radii [pm].
# Source: Kittel, "Introduction to Solid State Physics", 8th ed., Table 1.
_ATOMIC_RADII_PM = {
    "Fe": 126.0,
    "Co": 125.0,
    "Ni": 124.0,
    "Cu": 128.0,
    "Zn": 134.0,
    "Ga": 135.0,
    "In": 167.0,
    "Mo": 139.0,
    "W":  139.0,
    "Sn": 151.0,
}

# Valence electron count (VEC) per element.
# Convention (Guo & Liu, Prog. Nat. Sci. 2011): transition metals use the full
# group number (s + d electrons); p-block metals use s + p valence electrons.
# VEC of a composition is the atomic-fraction-weighted mean of these values.
# Add any missing element here to enable the VEC check for that composition.
_VALENCE_ELECTRONS = {
    "Fe": 8.0,
    "Co": 9.0,
    "Ni": 10.0,
    "Cu": 11.0,
    "Zn": 12.0,
    "Ga": 3.0,
    "In": 3.0,
    "Mo": 6.0,
    "W":  6.0,
    "Sn": 4.0,
}

# Binary mixing enthalpies ΔH^{A-B}_mix [kJ mol⁻¹].
# Values are estimated from Miedema's model or taken from the literature.
# Replace each None with the correct numerical value before use; pairs whose
# value is None cause the mixing-enthalpy check to be skipped (with a warning).
#
# Keys are alphabetically-sorted (A, B) tuples so both orderings resolve to the
# same entry automatically.
_BINARY_HMIX = {
    # Co — Fe
    ("Co", "Fe"):   -0.06904243624999928,   # ΔH_mix(Co–Fe)  [kJ/mol]
    # Co — Ni
    ("Co", "Ni"):   -0.02639796000000061,   # ΔH_mix(Co–Ni)  [kJ/mol]
    # Co — Cu
    ("Co", "Cu"):   0.05296564250000024,    # ΔH_mix(Co–Cu)  [kJ/mol]
    # Co — Zn
    ("Co", "Zn"):   -0.06440434124999994,   # ΔH_mix(Co–Zn)  [kJ/mol]
    # Co — Ga
    ("Co", "Ga"):   -0.2825187887500009,    # ΔH_mix(Co–Ga)  [kJ/mol]
    # Co — In
    ("Co", "In"):   -0.03301798562499947,   # ΔH_mix(Co–In)  [kJ/mol]
    # Co — Mo
    ("Co", "Mo"):   -0.044348996250000994,  # ΔH_mix(Co–Mo)  [kJ/mol]
    # Co — W
    ("Co", "W"):    -0.08396261250000059,   # ΔH_mix(Co–W)   [kJ/mol]
    # Co — Sn
    ("Co", "Sn"):   -0.14228259083333347,   # ΔH_mix(Co–Sn)  [kJ/mol]
    # Fe — Ni
    ("Fe", "Ni"):   -0.09145708749999937,   # ΔH_mix(Fe–Ni)  [kJ/mol]
    # Fe — Cu
    ("Cu", "Fe"):   0.06767582333333325,    # ΔH_mix(Cu–Fe)  [kJ/mol]
    # Fe — Zn
    ("Fe", "Zn"):   -0.03887403857142857,   # ΔH_mix(Fe–Zn)  [kJ/mol]
    # Fe — Ga
    ("Fe", "Ga"):   -0.23574756187500068,   # ΔH_mix(Fe–Ga)  [kJ/mol]
    # Fe — In
    ("Fe", "In"):   0.08910415999999977,   # ΔH_mix(Fe–In)  [kJ/mol]
    # Fe — Mo
    ("Fe", "Mo"):   -0.0030923849999998274,   # ΔH_mix(Fe–Mo)  [kJ/mol]
    # Fe — W
    ("Fe", "W"):    -0.023666749166667483,   # ΔH_mix(Fe–W)   [kJ/mol]
    # Fe — Sn
    ("Fe", "Sn"):   -0.0393573491666667,   # ΔH_mix(Fe–Sn)  [kJ/mol]
    # Ni — Cu
    ("Cu", "Ni"):   -0.0018628412500003577,   # ΔH_mix(Cu–Ni)  [kJ/mol]
    # Ni — Zn
    ("Ni", "Zn"):   -0.2547110099999994,   # ΔH_mix(Ni–Zn)  [kJ/mol]
    # Ni — Ga
    ("Ga", "Ni"):   -0.4091984924999991,   # ΔH_mix(Ga–Ni)  [kJ/mol]
    # Ni — In
    ("In", "Ni"):   -0.1919506283333329,   # ΔH_mix(In–Ni)  [kJ/mol]
    # Ni — Mo
    ("Mo", "Ni"):   -0.0923721724999993,   # ΔH_mix(Mo–Ni)  [kJ/mol]
    # Ni — W
    ("Ni", "W"):    -0.10733122500000025,   # ΔH_mix(Ni–W)   [kJ/mol]
    # Ni — Sn
    ("Ni", "Sn"):   -0.2849230889999994,   # ΔH_mix(Ni–Sn)  [kJ/mol]
    # Cu — Zn
    ("Cu", "Zn"):   -0.11121758769230758,   # ΔH_mix(Cu–Zn)  [kJ/mol]
    # Cu — Ga
    ("Cu", "Ga"):   -0.11059130346153814,   # ΔH_mix(Cu–Ga)  [kJ/mol]
    # Cu — In
    ("Cu", "In"):   -0.015922263333333575,   # ΔH_mix(Cu–In)  [kJ/mol]
    # Cu — Mo
    ("Cu", "Mo"):   0.07922715333333367,   # ΔH_mix(Cu–Mo)  [kJ/mol]
    # Cu — W
    ("Cu", "W"):    0.1272672366666671,   # ΔH_mix(Cu–W)   [kJ/mol]
    # Cu — Sn
    ("Cu", "Sn"):   -0.05965081250000015,   # ΔH_mix(Cu–Sn)  [kJ/mol]
    # Zn — Ga
    ("Ga", "Zn"):   0.014705678750000098,   # ΔH_mix(Ga–Zn)  [kJ/mol]
    # Zn — In
    ("In", "Zn"):   0.01643934833333353,   # ΔH_mix(In–Zn)  [kJ/mol]
    # Zn — Mo
    ("Mo", "Zn"):   -0.04813484562499992,   # ΔH_mix(Mo–Zn)  [kJ/mol]
    # Zn — W
    ("W",  "Zn"):   0.05189458250000012,   # ΔH_mix(W–Zn)   [kJ/mol]
    # Zn — Sn
    ("Sn", "Zn"):   0.02516417500000001,   # ΔH_mix(Sn–Zn)  [kJ/mol]
    # Ga — In
    ("Ga", "In"):   0.02130772708333299,   # ΔH_mix(Ga–In)  [kJ/mol]
    # Ga — Mo
    ("Ga", "Mo"):   -0.17489470199999957,   # ΔH_mix(Ga–Mo)  [kJ/mol]
    # Ga — W
    ("Ga", "W"):    -0.08640789800000022,   # ΔH_mix(Ga–W)   [kJ/mol]
    # Ga — Sn
    ("Ga", "Sn"):   0.035810431666666726,   # ΔH_mix(Ga–Sn)  [kJ/mol]
    # In — Mo
    ("In", "Mo"):   0.022434562000000113,   # ΔH_mix(In–Mo)  [kJ/mol]
    # In — W
    ("In", "W"):    0.1469121500000007,   # ΔH_mix(In–W)   [kJ/mol]
    # In — Sn
    ("In", "Sn"):   0.001599915625000392,   # ΔH_mix(In–Sn)  [kJ/mol]
    # Mo — W
    ("Mo", "W"):    -0.017974743333333265,   # ΔH_mix(Mo–W)   [kJ/mol]
    # Mo — Sn
    ("Mo", "Sn"):   -0.03018114499999945,   # ΔH_mix(Mo–Sn)  [kJ/mol]
    # W  — Sn
    ("Sn", "W"):    0.16359640333333422,   # ΔH_mix(Sn–W)   [kJ/mol]
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


def surface_id_from_composition(composition):
    # Preserve user-specified element order; it usually follows the alloy name.
    return "".join(f"{el}{count}" for el, count in composition)


def composition_from_cif(path):
    """Return deterministic element atom counts from a CIF."""
    atoms = read_cif(path)
    counts = Counter(atoms.get_chemical_symbols())
    if not counts:
        raise SystemExit(f"[error] CIF contains no atoms: {path}")
    return [(element, counts[element]) for element in sorted(counts)]


def compositions_match(left, right):
    return dict(left) == dict(right)


def dataset_paths(root, surface_id=None):
    root = Path(root)
    paths = {
        "root": root,
        "db": root / "index.sqlite",
        "manifest": root / "dataset_manifest.json",
    }
    if surface_id:
        surface = root / surface_id
        paths.update({
            "surface": surface,
            "structures": surface / "structures",
            "metadata": surface / "metadata",
            "openmx_slab": surface / "openmx_slab",
            "adsorbates": surface / "adsorbates",
            "surface_manifest": surface / "manifest.json",
        })
    return paths


def connect_db(root):
    return sqlite3.connect(dataset_paths(root)["db"])


def init_db(root):
    paths = dataset_paths(root)
    paths["root"].mkdir(parents=True, exist_ok=True)
    with sqlite3.connect(paths["db"]) as con:
        old_version = con.execute(
            "SELECT value FROM dataset_info WHERE key='schema_version'"
        ).fetchone() if con.execute(
            "SELECT 1 FROM sqlite_master WHERE type='table' AND name='dataset_info'"
        ).fetchone() else None
        if old_version and int(old_version[0]) != SCHEMA_VERSION:
            raise SystemExit(
                f"[error] Dataset schema v{old_version[0]} is incompatible with v{SCHEMA_VERSION}. "
                "Use a new dataset root or migrate the existing dataset explicitly."
            )
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
                path TEXT NOT NULL,
                status TEXT NOT NULL DEFAULT 'created',
                initial_cif TEXT,
                relaxed_cif TEXT,
                total_energy_eV REAL,
                created_at TEXT NOT NULL,
                updated_at TEXT NOT NULL
            );

            CREATE TABLE IF NOT EXISTS artifacts (
                artifact_id INTEGER PRIMARY KEY AUTOINCREMENT,
                surface_id TEXT NOT NULL,
                kind TEXT NOT NULL,
                path TEXT NOT NULL,
                sha256 TEXT,
                created_at TEXT NOT NULL,
                FOREIGN KEY (surface_id) REFERENCES surfaces(surface_id)
            );

            CREATE TABLE IF NOT EXISTS top_atoms (
                surface_id TEXT NOT NULL,
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
                PRIMARY KEY (surface_id, atom_id),
                FOREIGN KEY (surface_id) REFERENCES surfaces(surface_id)
            );

            CREATE TABLE IF NOT EXISTS fcc_sites (
                surface_id TEXT NOT NULL,
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
                PRIMARY KEY (surface_id, site_id),
                FOREIGN KEY (surface_id) REFERENCES surfaces(surface_id)
            );

            CREATE TABLE IF NOT EXISTS adsorbate_configs (
                config_id TEXT PRIMARY KEY,
                surface_id TEXT NOT NULL,
                site_id TEXT NOT NULL,
                adsorbate TEXT NOT NULL,
                path TEXT NOT NULL,
                initial_cif TEXT,
                relaxed_cif TEXT,
                adsorption_energy_eV REAL,
                energy_status TEXT NOT NULL DEFAULT 'empty',
                updated_at TEXT NOT NULL,
                FOREIGN KEY (surface_id, site_id)
                    REFERENCES fcc_sites(surface_id, site_id)
            );

            CREATE TABLE IF NOT EXISTS hamiltonian_exports (
                export_id INTEGER PRIMARY KEY AUTOINCREMENT,
                surface_id TEXT NOT NULL,
                scfout_path TEXT NOT NULL,
                output_npz TEXT NOT NULL,
                n_surface_atoms INTEGER NOT NULL,
                spin_channels INTEGER NOT NULL,
                created_at TEXT NOT NULL,
                FOREIGN KEY (surface_id) REFERENCES surfaces(surface_id)
            );
            """
        )
        # Additive migration: add total_energy_eV to existing datasets without
        # bumping SCHEMA_VERSION (nullable column, fully backwards-compatible).
        try:
            con.execute("ALTER TABLE surfaces ADD COLUMN total_energy_eV REAL")
        except sqlite3.OperationalError:
            pass  # column already present
        con.execute(
            "INSERT OR REPLACE INTO dataset_info(key, value) VALUES (?, ?)",
            ("schema_version", str(SCHEMA_VERSION)),
        )
        con.commit()

    write_json(paths["manifest"], {
        "schema_version": SCHEMA_VERSION,
        "created_or_updated_at": utc_now(),
        "layout": "<surface_id>",
        "surface_id_rule": "concatenated_element_atom_count_pairs_only",
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


def surface_exists(root, surface_id):
    """Return True if the surface is registered in the DB or present on disk."""
    paths = dataset_paths(root, surface_id)
    if paths["surface"].exists():
        return True
    db_path = dataset_paths(root)["db"]
    if db_path.is_file():
        with sqlite3.connect(db_path) as con:
            has_table = con.execute(
                "SELECT 1 FROM sqlite_master WHERE type='table' AND name='surfaces'"
            ).fetchone()
            if has_table and con.execute(
                "SELECT 1 FROM surfaces WHERE surface_id=?", (surface_id,)
            ).fetchone():
                return True
    return False


def composition_from_counts(elements, counts):
    """Build a sorted (element, count) composition, dropping zero-count elements.

    Sorting alphabetically and dropping zeros makes the surface_id match the one
    inferred later from a CIF via composition_from_cif, so create-sample and
    record-relaxed agree on the identity of a composition.
    """
    merged = Counter()
    for element, count in zip(elements, counts):
        merged[element] += int(count)
    return [(element, merged[element]) for element in sorted(merged)
            if merged[element] > 0]


def _check_hea_criteria(elements, fractions):
    """Screen a composition against three HEA formation criteria.

    Args:
        elements: list of element symbols
        fractions: normalized atomic fractions (mole fractions)

    Returns:
        (passed: bool, report: dict)
            passed is True when all computable checks pass.
            report contains per-criterion values and pass/fail status.
    """
    xs = np.array(fractions, dtype=float)
    xs = xs / xs.sum()  # renormalize for safety

    # -----------------------------------------------------------------------
    # 1. Configurational entropy:  ΔS_mix = -R Σ x_i ln(x_i)
    #    Criterion: ΔS_mix > 1.5 R
    # -----------------------------------------------------------------------
    nonzero = xs > 1e-12
    s_mix = -_R * np.sum(xs[nonzero] * np.log(xs[nonzero]))
    entropy_ok = bool(s_mix > _HEA_ENTROPY_MIN)

    # -----------------------------------------------------------------------
    # 2. Atomic size mismatch:  δ = 100 × sqrt(Σ x_i (1 - r_i / r̄)²)
    #    Criterion: δ ≤ 6.6 %
    # -----------------------------------------------------------------------
    radii = np.array([_ATOMIC_RADII_PM.get(el) for el in elements])
    missing_radii = [el for el, r in zip(elements, radii) if r is None]
    if missing_radii:
        delta = None
        size_ok = None
    else:
        radii = radii.astype(float)
        r_bar = np.sum(xs * radii)
        delta = 100.0 * np.sqrt(np.sum(xs * ((1.0 - radii / r_bar) ** 2)))
        size_ok = bool(delta <= _HEA_SIZE_DELTA_MAX)

    # -----------------------------------------------------------------------
    # 3. Mixing enthalpy:  ΔH_mix ≈ Σ_{i<j} 4·ΔH^{ij}_binary · x_i · x_j
    #    Criterion: -15 kJ/mol ≤ ΔH_mix ≤ +5 kJ/mol
    # -----------------------------------------------------------------------
    h_mix = 0.0
    missing_pairs = []
    can_compute_h = True
    for i, (el_i, x_i) in enumerate(zip(elements, xs)):
        for j in range(i + 1, len(elements)):
            el_j = elements[j]
            x_j = xs[j]
            key = tuple(sorted([el_i, el_j]))
            val = _BINARY_HMIX.get(key)
            if val is None:
                missing_pairs.append(f"{el_i}-{el_j}")
                can_compute_h = False
            else:
                h_mix += 4.0 * val * x_i * x_j

    if not can_compute_h:
        h_mix = None
        hmix_ok = None
    else:
        hmix_ok = bool(_HEA_HMIX_MIN <= h_mix <= _HEA_HMIX_MAX)

    # -----------------------------------------------------------------------
    # 4. Valence electron concentration:  VEC = Σ x_i · VEC_i
    #    Criterion: VEC ≥ 8
    # -----------------------------------------------------------------------
    vec_values = np.array([_VALENCE_ELECTRONS.get(el) for el in elements])
    missing_vec = [el for el, v in zip(elements, vec_values) if v is None]
    if missing_vec:
        vec = None
        vec_ok = None
    else:
        vec = float(np.sum(xs * vec_values.astype(float)))
        vec_ok = bool(vec >= _HEA_VEC_MIN)

    # All computable checks must pass
    passed = all(
        c is True
        for c in [entropy_ok, size_ok, hmix_ok, vec_ok]
        if c is not None
    )

    report = {
        "s_mix_over_R": float(s_mix / _R),
        "entropy_ok": entropy_ok,
        "delta_pct": float(delta) if delta is not None else None,
        "size_ok": size_ok,
        "h_mix_kJ_mol": float(h_mix) if h_mix is not None else None,
        "hmix_ok": hmix_ok,
        "vec": vec,
        "vec_ok": vec_ok,
        "missing_radii": missing_radii,
        "missing_hmix_pairs": missing_pairs,
        "missing_vec": missing_vec,
    }
    return passed, report


def create_sample(args):
    """Generate an SQS slab from elements (+optional ratios) and register it.

    When ratios are not given they are drawn at random; if the resulting
    composition already exists, new ratios are drawn until a novel composition
    is found (bounded by --max-attempts). Fixed user ratios that collide are a
    hard error because they cannot be regenerated.
    """
    try:
        import build_hea_surface as builder
    except ImportError as exc:
        raise SystemExit(f"[error] Could not import build_hea_surface.py: {exc}")

    elements = args.elements
    rng = np.random.default_rng(args.seed)
    size = tuple(args.size)
    random_ratios = args.ratios is None

    # Build the template once to learn the site count; substitution/SQS reuse it.
    a = args.lattice_constant
    template = None
    n_sites = None
    composition = None
    counts = None

    hea_report = None
    for attempt in range(1, args.max_attempts + 1):
        fractions = builder.normalize_ratios(
            elements, args.ratios, rng, verbose=(attempt == 1 and random_ratios)
        )
        if template is None:
            lattice = a or builder.estimate_lattice_constant(elements, fractions)
            template = builder.build_template_slab(
                elements[0], size, lattice, args.vacuum
            )
            n_sites = len(template)

        # --- HEA thermodynamic screening ---
        _, hea_report = _check_hea_criteria(elements, fractions)

        # Entropy, size mismatch, mixing enthalpy and VEC all depend on the
        # ratio, so we retry when ratios are random, or abort immediately for
        # fixed ratios. A criterion whose inputs are missing evaluates to None
        # and is skipped (treated as non-blocking), mirroring the mixing-enthalpy
        # behaviour for incomplete binary-pair data.
        entropy_ok = hea_report["entropy_ok"]
        size_ok = hea_report["size_ok"]  # None means skipped (missing radii)
        hmix_ok = hea_report["hmix_ok"]  # None means skipped (missing pairs)
        vec_ok = hea_report["vec_ok"]    # None means skipped (missing VEC data)
        ratio_criteria_ok = (
            entropy_ok
            and (size_ok is not False)
            and (hmix_ok is not False)
            and (vec_ok is not False)
        )
        if not ratio_criteria_ok:
            if not random_ratios:
                raise SystemExit(
                    "[error] Specified composition does not satisfy HEA criteria:\n"
                    f"  ΔS_mix/R = {hea_report['s_mix_over_R']:.3f}  "
                    f"(need > {_HEA_ENTROPY_MIN / _R:.1f})  "
                    f"{'OK' if entropy_ok else 'FAIL'}\n"
                    f"  δ        = "
                    + (f"{hea_report['delta_pct']:.2f}%" if hea_report["delta_pct"] is not None else "N/A")
                    + f"  (need ≤ {_HEA_SIZE_DELTA_MAX}%)  "
                    + ("OK" if size_ok is True else ("FAIL" if size_ok is False else "SKIP"))
                    + "\n"
                    f"  ΔH_mix   = "
                    + (f"{hea_report['h_mix_kJ_mol']:.2f} kJ/mol" if hea_report["h_mix_kJ_mol"] is not None else "N/A")
                    + f"  (need [{_HEA_HMIX_MIN}, {_HEA_HMIX_MAX}] kJ/mol)  "
                    + ("OK" if hmix_ok is True else ("FAIL" if hmix_ok is False else "SKIP"))
                    + "\n"
                    f"  VEC      = "
                    + (f"{hea_report['vec']:.2f}" if hea_report["vec"] is not None else "N/A")
                    + f"  (need ≥ {_HEA_VEC_MIN:.1f})  "
                    + ("OK" if vec_ok is True else ("FAIL" if vec_ok is False else "SKIP"))
                )
            if attempt == args.max_attempts:
                raise SystemExit(
                    f"[error] Could not find a composition satisfying the entropy, "
                    f"size-mismatch, mixing-enthalpy and VEC HEA criteria after "
                    f"{args.max_attempts} attempts."
                )
            continue  # draw new random fractions
        # -----------------------------------

        candidate_counts = builder.largest_remainder_counts(fractions, n_sites)
        candidate = composition_from_counts(elements, candidate_counts)
        surface_id = surface_id_from_composition(candidate)
        if not surface_exists(args.root, surface_id):
            composition = candidate
            counts = candidate_counts
            break
        if not random_ratios:
            raise SystemExit(
                f"[error] Composition already exists in dataset: {surface_id}"
            )
        if attempt == args.max_attempts:
            raise SystemExit(
                f"[error] Could not find a novel composition after "
                f"{args.max_attempts} attempts; the composition space for these "
                f"elements and {n_sites} sites may be exhausted."
            )

    surface_id = surface_id_from_composition(composition)
    paths = dataset_paths(args.root, surface_id)

    print("=" * 60)
    print(f"Template      : FCC(111), size={size}, sites={n_sites}")
    print(f"Surface ID    : {surface_id}")
    print("Composition   :")
    for element, count in composition:
        print(f"  {element:>3s} : {count:>3d} atoms ({count / n_sites * 100:6.2f}%)")
    if hea_report is not None:
        print("HEA criteria  :")
        print(
            f"  ΔS_mix/R = {hea_report['s_mix_over_R']:.3f}"
            f"  (threshold > {_HEA_ENTROPY_MIN / _R:.1f})  "
            f"{'[OK]' if hea_report['entropy_ok'] else '[FAIL]'}"
        )
        if hea_report["delta_pct"] is not None:
            print(
                f"  δ        = {hea_report['delta_pct']:.2f}%"
                f"  (threshold ≤ {_HEA_SIZE_DELTA_MAX}%)  "
                f"{'[OK]' if hea_report['size_ok'] else '[FAIL]'}"
            )
        else:
            print(f"  δ        = N/A  (missing radii: {', '.join(hea_report['missing_radii'])})")
        if hea_report["h_mix_kJ_mol"] is not None:
            print(
                f"  ΔH_mix   = {hea_report['h_mix_kJ_mol']:.2f} kJ/mol"
                f"  (threshold [{_HEA_HMIX_MIN}, {_HEA_HMIX_MAX}] kJ/mol)  "
                f"{'[OK]' if hea_report['hmix_ok'] else '[FAIL]'}"
            )
        else:
            print(f"  ΔH_mix   = N/A  (missing pairs: {', '.join(hea_report['missing_hmix_pairs'])})")
        if hea_report["vec"] is not None:
            print(
                f"  VEC      = {hea_report['vec']:.2f}"
                f"  (threshold ≥ {_HEA_VEC_MIN:.1f})  "
                f"{'[OK]' if hea_report['vec_ok'] else '[FAIL]'}"
            )
        else:
            print(f"  VEC      = N/A  (missing VEC data: {', '.join(hea_report['missing_vec'])})")
    print("=" * 60)

    # Random substitution -> initial structure, then optional SQS refinement.
    random_slab = builder.assign_random_substitution(
        template, elements, counts, rng
    )
    final = None
    if not args.no_sqs:
        final = builder.run_sqs(
            template, elements, counts, args.cutoffs, args.n_steps, args.seed
        )
    if final is None:
        final = random_slab
        method = "random"
    else:
        method = "sqs"
        print("[info] SQS optimization finished.")
    final = builder.sort_atoms_by_element(final)

    init_db(args.root)
    for key in ("structures", "metadata", "openmx_slab", "adsorbates"):
        paths[key].mkdir(parents=True, exist_ok=True)

    initial_cif = paths["structures"] / "00_initial_sqs.cif"
    from ase.io import write as ase_write
    ase_write(str(initial_cif), final)

    manifest = {
        "surface_id": surface_id,
        "composition": [{"element": el, "atom_count": count}
                        for el, count in composition],
        "generation": {
            "method": method,
            "size": list(size),
            "n_sites": n_sites,
            "vacuum": args.vacuum,
            "random_ratios": random_ratios,
            "seed": args.seed,
        },
        "created_at": utc_now(),
        "files": {
            "initial_cif": relpath_or_none(initial_cif, paths["surface"]),
            "relaxed_cif": None,
            "top_atoms": "metadata/top_atoms.jsonl",
            "fcc_sites": "metadata/fcc_sites.jsonl",
            "atom_grid": "metadata/atom_grid.npy",
            "site_grid": "metadata/site_grid.npy",
        },
    }
    write_json(paths["surface_manifest"], manifest)

    with connect_db(args.root) as con:
        con.execute(
            """INSERT INTO surfaces(
                surface_id, composition_json, path, status, initial_cif,
                relaxed_cif, created_at, updated_at
            ) VALUES (?, ?, ?, ?, ?, ?, ?, ?)""",
            (surface_id, json.dumps(manifest["composition"]), str(paths["surface"]),
             "created", relpath_or_none(initial_cif, paths["surface"]),
             None, utc_now(), utc_now()),
        )
        con.execute(
            """
            INSERT INTO artifacts(surface_id, kind, path, sha256, created_at)
            VALUES (?, ?, ?, ?, ?)
            """,
            (surface_id, "initial_cif", str(initial_cif),
             sha256_file(initial_cif), utc_now()),
        )
        con.commit()

    print(f"[done] Created composition {surface_id} ({method}) -> {initial_cif}")


def record_relaxed(args):
    """Register an externally relaxed slab against an existing composition."""
    paths = dataset_paths(args.root, args.surface_id)
    initial_cif = paths["structures"] / "00_initial_sqs.cif"
    if not surface_exists(args.root, args.surface_id) or not initial_cif.is_file():
        raise SystemExit(
            f"[error] Surface not registered; run create-sample first: "
            f"{args.surface_id}"
        )

    relaxed_source = Path(args.relaxed_cif)
    if not relaxed_source.is_file():
        raise SystemExit(f"[error] Relaxed CIF not found: {relaxed_source}")

    initial_composition = composition_from_cif(initial_cif)
    relaxed_composition = composition_from_cif(relaxed_source)
    if not compositions_match(initial_composition, relaxed_composition):
        raise SystemExit(
            "[error] Relaxed slab composition differs from the registered "
            f"initial structure: {surface_id_from_composition(relaxed_composition)}"
        )

    # index-surface associates relaxed coordinates by ASE atom index, so the
    # relaxed CIF must preserve both atom count and per-atom element order.
    initial_atoms = read_cif(initial_cif)
    relaxed_atoms = read_cif(relaxed_source)
    if len(initial_atoms) != len(relaxed_atoms):
        raise SystemExit(
            f"[error] Atom count changed during relaxation: "
            f"{len(initial_atoms)} -> {len(relaxed_atoms)}."
        )
    if initial_atoms.get_chemical_symbols() != relaxed_atoms.get_chemical_symbols():
        raise SystemExit(
            "[error] Relaxed slab CIF does not preserve the initial atom order."
        )

    relaxed_target = paths["structures"] / "01_relaxed_slab.cif"
    if relaxed_source.resolve() != relaxed_target.resolve():
        shutil.copy2(relaxed_source, relaxed_target)

    if paths["surface_manifest"].is_file():
        with open(paths["surface_manifest"], "r", encoding="utf-8") as f:
            manifest = json.load(f)
        manifest.setdefault("files", {})["relaxed_cif"] = relpath_or_none(
            relaxed_target, paths["surface"]
        )
        if args.energy is not None:
            manifest["total_energy_eV"] = args.energy
        write_json(paths["surface_manifest"], manifest)

    with connect_db(args.root) as con:
        con.execute(
            "UPDATE surfaces SET relaxed_cif=?, status=?, total_energy_eV=?, updated_at=? "
            "WHERE surface_id=?",
            (relpath_or_none(relaxed_target, paths["surface"]), "slab_relaxed",
             args.energy, utc_now(), args.surface_id),
        )
        con.execute(
            "DELETE FROM artifacts WHERE surface_id=? AND kind='relaxed_cif'",
            (args.surface_id,),
        )
        con.execute(
            """
            INSERT INTO artifacts(surface_id, kind, path, sha256, created_at)
            VALUES (?, ?, ?, ?, ?)
            """,
            (args.surface_id, "relaxed_cif", str(relaxed_target),
             sha256_file(relaxed_target), utc_now()),
        )
        con.commit()

    print(f"[done] Recorded relaxed slab for {args.surface_id} -> {relaxed_target}")


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
    paths = dataset_paths(args.root, args.surface_id)
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
            "DELETE FROM top_atoms WHERE surface_id=?", (args.surface_id,),
        )
        for item in rows:
            con.execute(
                """
                INSERT INTO top_atoms(
                    surface_id, atom_id, row, col, element,
                    initial_ase_index, relaxed_ase_index, openmx_atom_index,
                    initial_x, initial_y, initial_z, relaxed_x, relaxed_y, relaxed_z
                ) VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)
                """,
                (
                    item["surface_id"], item["atom_id"],
                    item["row"], item["col"], item["element"],
                    item["initial_ase_index"], item["relaxed_ase_index"],
                    item["openmx_atom_index"], item["initial_x"], item["initial_y"],
                    item["initial_z"], item["relaxed_x"], item["relaxed_y"],
                    item["relaxed_z"],
                ),
            )
        con.execute(
            "UPDATE surfaces SET status=?, updated_at=? WHERE surface_id=?",
            ("surface_indexed", utc_now(), args.surface_id),
        )
        con.commit()

    print(f"[done] Wrote {top_atoms_path} ({len(rows)} top atoms)")


def detect_sites(args):
    try:
        import add_fcc_adsorbate as ads
    except ImportError as exc:
        raise SystemExit(f"[error] Could not import add_fcc_adsorbate.py: {exc}")

    paths = dataset_paths(args.root, args.surface_id)
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
            "DELETE FROM fcc_sites WHERE surface_id=?", (args.surface_id,),
        )
        for item in rows:
            con.execute(
                """
                INSERT INTO fcc_sites(
                    surface_id, site_id, site_index, row, col,
                    site_type, frac_x, frac_y, plane_x, plane_y, plane_z,
                    top_atom_ids_json, adsorption_energy_eV, energy_status
                ) VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)
                """,
                (
                    item["surface_id"], item["site_id"],
                    item["site_index"], item["row"], item["col"],
                    item["site_type"], item["frac_x"], item["frac_y"],
                    item["plane_x"], item["plane_y"], item["plane_z"],
                    json.dumps(item["top_atom_ids"]),
                    item["adsorption_energy_eV"], item["energy_status"],
                ),
            )
        con.execute(
            "UPDATE surfaces SET status=?, updated_at=? WHERE surface_id=?",
            ("sites_detected", utc_now(), args.surface_id),
        )
        con.commit()

    print(f"[done] Wrote {fcc_path} ({len(rows)} FCC sites)")


def create_adsorbate_records(args):
    try:
        import add_fcc_adsorbate as ads
    except ImportError as exc:
        raise SystemExit(f"[error] Could not import add_fcc_adsorbate.py: {exc}")

    paths = dataset_paths(args.root, args.surface_id)
    relaxed_cif = paths["structures"] / "01_relaxed_slab.cif"
    if not relaxed_cif.is_file():
        raise SystemExit(
            f"[error] Registered relaxed slab CIF not found: {relaxed_cif}"
        )
    slab = ads.read_cif(str(relaxed_cif))
    slab.set_pbc(True)
    normal = ads.normal_from_cell(slab.cell, args.side)

    fcc_rows = read_jsonl(paths["metadata"] / "fcc_sites.jsonl")
    if not fcc_rows:
        raise SystemExit("[error] No FCC site records found; run detect-sites first.")
    ads_root = paths["adsorbates"] / args.adsorbate
    ads_root.mkdir(parents=True, exist_ok=True)

    rebuild = getattr(args, "rebuild_initial", False)
    created = 0
    preserved = 0
    with connect_db(args.root) as con:
        for site in fcc_rows:
            site_id = site["site_id"]
            config_id = f"{args.surface_id}_{args.adsorbate}_{site_id}"
            config_path = ads_root / site_id
            initial_adsorbate_cif = config_path / "00_initial_adsorbate.cif"

            # Non-destructive resume: never clobber a config that already exists,
            # since it may carry a computed adsorption energy. Only an explicit
            # --rebuild-initial regenerates the structure and resets the status.
            already = con.execute(
                "SELECT 1 FROM adsorbate_configs WHERE config_id=?", (config_id,)
            ).fetchone() is not None
            if not rebuild and already and initial_adsorbate_cif.is_file():
                preserved += 1
                continue

            config_path.mkdir(parents=True, exist_ok=True)
            site_geometry = {
                "plane_pos": np.array(
                    [site["plane_x"], site["plane_y"], site["plane_z"]],
                    dtype=float,
                )
            }
            adsorbed = ads.add_adsorbate(
                slab, site_geometry, args.adsorbate, normal, args.height, args.nh
            )
            ads.write_cif(str(initial_adsorbate_cif), adsorbed)
            record = {
                "surface_id": args.surface_id,
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
                    config_id, surface_id, site_id, adsorbate, path,
                    initial_cif, relaxed_cif, adsorption_energy_eV, energy_status,
                    updated_at
                ) VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?)
                """,
                (
                    config_id, args.surface_id, site_id, args.adsorbate, str(config_path),
                    "00_initial_adsorbate.cif", "01_relaxed_adsorbate.cif",
                    None, "empty", utc_now(),
                ),
            )
            created += 1
        con.commit()
    print(
        f"[done] Adsorbate records under {ads_root} "
        f"(created/reset {created}, preserved {preserved})"
    )


def validate_relaxed_adsorbate(initial_path, relaxed_path, clean_slab_path,
                                max_displacement):
    """Validate atom identity and adsorbate displacement after relaxation."""
    from ase.geometry import find_mic

    initial = read_cif(initial_path)
    relaxed = read_cif(relaxed_path)
    clean_slab = read_cif(clean_slab_path)
    initial_symbols = initial.get_chemical_symbols()
    relaxed_symbols = relaxed.get_chemical_symbols()
    if len(initial) != len(relaxed):
        raise SystemExit(
            f"[error] Atom count changed during relaxation: "
            f"{len(initial)} -> {len(relaxed)}."
        )
    if initial_symbols != relaxed_symbols:
        raise SystemExit(
            "[error] Relaxed adsorbate CIF does not preserve the initial atom order."
        )

    n_slab = len(clean_slab)
    if n_slab >= len(initial):
        raise SystemExit("[error] Initial adsorbate CIF contains no adsorbate atoms.")
    if initial_symbols[:n_slab] != clean_slab.get_chemical_symbols():
        raise SystemExit(
            "[error] Initial adsorbate CIF substrate does not match the registered slab."
        )

    displacement_vectors, _ = find_mic(
        relaxed.positions - initial.positions,
        cell=initial.cell,
        pbc=initial.pbc,
    )
    # Remove a rigid translation of the whole relaxed structure using the
    # component-wise median displacement of substrate atoms.
    substrate_shift = np.median(displacement_vectors[:n_slab], axis=0)
    adsorbate_vectors = displacement_vectors[n_slab:] - substrate_shift
    adsorbate_displacements = np.linalg.norm(adsorbate_vectors, axis=1)
    maximum = float(np.max(adsorbate_displacements))
    if maximum > max_displacement:
        raise SystemExit(
            f"[error] Adsorbate moved too far during relaxation: {maximum:.4f} A "
            f"> allowed {max_displacement:.4f} A."
        )
    return {
        "max_adsorbate_displacement_A": maximum,
        "adsorbate_displacements_A": [
            float(value) for value in adsorbate_displacements
        ],
        "max_allowed_adsorbate_displacement_A": float(max_displacement),
    }


def record_energy(args):
    if args.max_adsorbate_displacement < 0:
        raise SystemExit("[error] --max-adsorbate-displacement must be non-negative.")
    if args.site_id < 1:
        raise SystemExit("[error] --site-id must be a positive integer.")
    site_id = f"site_{args.site_id:04d}"
    paths = dataset_paths(args.root, args.surface_id)
    config_path = paths["adsorbates"] / args.adsorbate / site_id
    energy_path = config_path / "adsorption_energy.json"
    if not energy_path.is_file():
        raise SystemExit(f"[error] Missing energy record: {energy_path}")
    initial_path = config_path / "00_initial_adsorbate.cif"
    clean_slab_path = paths["structures"] / "01_relaxed_slab.cif"
    relaxed_source = Path(args.relaxed_cif)
    validation = validate_relaxed_adsorbate(
        initial_path, relaxed_source, clean_slab_path,
        args.max_adsorbate_displacement,
    )
    relaxed_target = config_path / "01_relaxed_adsorbate.cif"
    if relaxed_source.resolve() != relaxed_target.resolve():
        shutil.copy2(relaxed_source, relaxed_target)

    with open(energy_path, "r", encoding="utf-8") as f:
        record = json.load(f)
    record["relaxed_cif"] = "01_relaxed_adsorbate.cif"
    record["adsorption_energy_eV"] = args.energy
    record["energy_status"] = args.status
    record["relaxation_validation"] = validation
    record["notes"] = args.notes or record.get("notes", "")
    write_json(energy_path, record)

    with connect_db(args.root) as con:
        con.execute(
            """
            UPDATE adsorbate_configs
            SET adsorption_energy_eV=?, energy_status=?, updated_at=?
            WHERE surface_id=? AND adsorbate=? AND site_id=?
            """,
            (
                args.energy, args.status, utc_now(), args.surface_id,
                args.adsorbate, site_id,
            ),
        )
        con.execute(
            """
            UPDATE fcc_sites
            SET adsorption_energy_eV=?, energy_status=?
            WHERE surface_id=? AND site_id=?
            """,
            (args.energy, args.status, args.surface_id, site_id),
        )
        con.commit()
    print(f"[done] Recorded {args.energy} eV for {args.adsorbate}/{site_id}")


def extract_hamiltonian(args):
    """Extract top-layer d-orbital Hamiltonian data for a registered surface."""
    try:
        import extract_openmx_hamiltonian as ex
    except ImportError as exc:
        raise SystemExit(
            f"[error] Could not import extract_openmx_hamiltonian.py: {exc}"
        )

    paths = dataset_paths(args.root, args.surface_id)
    top_atoms = paths["metadata"] / "top_atoms.jsonl"
    if not top_atoms.is_file():
        raise SystemExit(
            f"[error] Top-atom metadata not found: {top_atoms}; run index-surface first."
        )
    dat_path = Path(args.dat)
    scfout_path = Path(args.scfout)
    output = Path(
        args.output or paths["openmx_slab"] / "hamiltonian_d_surface.npz"
    )

    species_basis, atoms = ex.parse_dat(dat_path)
    basis, offsets = ex.build_basis(atoms, species_basis)
    scfout = ex.parse_scfout_binary(
        scfout_path, expected_total_orbitals=len(basis)
    )
    top_rows = read_jsonl(top_atoms)
    atom_ids, openmx_indices, d_lists, d_label_lists = ex.d_indices_for_surface(
        top_rows, basis, offsets
    )
    h_d, d_basis = ex.extract_blocks(scfout.hamiltonian, d_lists)

    output.parent.mkdir(parents=True, exist_ok=True)
    max_d = d_basis.shape[1]
    d_labels = np.full((len(d_label_lists), max_d), "", dtype="<U32")
    for index, labels in enumerate(d_label_lists):
        d_labels[index, :len(labels)] = labels
    np.savez_compressed(
        output,
        H_d=h_d,
        d_basis_indices=d_basis,
        d_labels=d_labels,
        surface_atom_ids=np.array(atom_ids, dtype="<U64"),
        openmx_atom_indices=np.array(openmx_indices, dtype=np.int64),
        spin_switch=np.array([scfout.spin_switch], dtype=np.int64),
        source_scfout=np.array([str(scfout_path)], dtype="<U1024"),
        source_dat=np.array([str(dat_path)], dtype="<U1024"),
    )
    basis_output = (
        Path(args.basis_output)
        if args.basis_output
        else output.with_suffix(output.suffix + ".basis.jsonl")
    )
    ex.write_basis_jsonl(basis_output, basis)

    with connect_db(args.root) as con:
        con.execute(
            """INSERT INTO hamiltonian_exports(
                surface_id, scfout_path, output_npz, n_surface_atoms,
                spin_channels, created_at
            ) VALUES (?, ?, ?, ?, ?, ?)""",
            (args.surface_id, str(scfout_path), str(output), len(atom_ids),
             h_d.shape[0], utc_now()),
        )
        con.commit()
    print(f"[done] Wrote {output}")
    print(f"[done] Wrote {basis_output}")
    print(f"[info] surface_atoms={len(atom_ids)} spin_channels={h_d.shape[0]}")


# ---------------------------------------------------------------------------
# Completeness check
# ---------------------------------------------------------------------------

_CHECK_TAG = {"ok": "[OK]     ", "missing": "[MISSING]", "partial": "[PARTIAL]"}


def _count_rows(con, table, surface_id):
    """Return the number of rows a surface owns in ``table`` (0 if absent)."""
    if con is None:
        return 0
    try:
        row = con.execute(
            f"SELECT COUNT(*) FROM {table} WHERE surface_id=?", (surface_id,)
        ).fetchone()
    except sqlite3.OperationalError:
        return 0  # table missing on an older dataset
    return int(row[0]) if row else 0


def _combine_check(label, db_present, disk_present, detail=None):
    """Fold DB-row presence and on-disk presence into one check outcome.

    ``ok`` when both agree the artifact exists, ``missing`` when both agree it
    does not, and ``partial`` when the index and the filesystem disagree (an
    orphan row or an unregistered file) so the inconsistency is surfaced rather
    than hidden.
    """
    if db_present and disk_present:
        return (label, "ok", detail)
    if not db_present and not disk_present:
        return (label, "missing", "not recorded")
    where = "index only" if db_present else "disk only"
    return (label, "partial", f"present in {where}")


def _fmt_id_list(ids, max_list):
    ids = list(ids)
    if max_list and len(ids) > max_list:
        return ", ".join(ids[:max_list]) + f", ... (+{len(ids) - max_list} more)"
    return ", ".join(ids)


def _check_adsorbate_coverage(con, paths, surface_id, adsorbate, site_ids):
    """Per-site coverage of one adsorbate: registered config, relaxed structure,
    and adsorption energy. Returns a dict summarising counts and missing sites."""
    total = len(site_ids)
    rows = con.execute(
        "SELECT site_id, relaxed_cif, adsorption_energy_eV "
        "FROM adsorbate_configs WHERE surface_id=? AND adsorbate=?",
        (surface_id, adsorbate),
    ).fetchall() if con is not None else []
    by_site = {r[0]: r for r in rows}

    missing_config, missing_struct, missing_energy = [], [], []
    for site_id in site_ids:
        row = by_site.get(site_id)
        if row is None:
            missing_config.append(site_id)
            missing_struct.append(site_id)
            missing_energy.append(site_id)
            continue
        relaxed_name = row[1] or "01_relaxed_adsorbate.cif"
        relaxed_cif = paths["adsorbates"] / adsorbate / site_id / relaxed_name
        if not relaxed_cif.is_file():
            missing_struct.append(site_id)
        if row[2] is None:
            missing_energy.append(site_id)

    complete = (
        total > 0
        and not missing_config and not missing_struct and not missing_energy
    )
    return {
        "adsorbate": adsorbate,
        "total_sites": total,
        "registered": total - len(missing_config),
        "structures": total - len(missing_struct),
        "energies": total - len(missing_energy),
        "missing_config": missing_config,
        "missing_structure": missing_struct,
        "missing_energy": missing_energy,
        "complete": complete,
    }


def check_surface(args):
    """Report completeness of one surface_id and optionally adsorbate coverage.

    Verifies the presence and DB/disk agreement of every workflow product:
    initial + relaxed structures, slab total energy, top-atom indexing, FCC
    site detection and the Hamiltonian export. When --adsorbates is given it
    additionally reports, per species, how many sites carry a registered
    config, a relaxed structure on disk and an adsorption energy in the index.

    Returns 0 when everything checked is complete, 1 otherwise, so callers can
    gate on the exit code.
    """
    root = args.root
    surface_id = args.surface_id
    if not surface_exists(root, surface_id):
        raise SystemExit(f"[error] Surface not found in dataset: {surface_id}")

    paths = dataset_paths(root, surface_id)
    db_path = dataset_paths(root)["db"]
    max_list = args.max_list if args.max_list and args.max_list > 0 else None

    con = sqlite3.connect(db_path) if db_path.is_file() else None
    try:
        surf_row = con.execute(
            "SELECT initial_cif, relaxed_cif, total_energy_eV, status "
            "FROM surfaces WHERE surface_id=?", (surface_id,)
        ).fetchone() if con is not None else None
        db_initial, db_relaxed, db_energy, status = (
            surf_row if surf_row else (None, None, None, None)
        )

        checks = []
        # 1. Initial HEA structure
        checks.append(_combine_check(
            "Initial structure registered",
            db_present=db_initial is not None,
            disk_present=(paths["structures"] / "00_initial_sqs.cif").is_file(),
        ))
        # 2. Relaxed HEA structure
        checks.append(_combine_check(
            "Relaxed structure registered",
            db_present=db_relaxed is not None,
            disk_present=(paths["structures"] / "01_relaxed_slab.cif").is_file(),
        ))
        # 3. Slab total energy (index-only field)
        checks.append((
            "Slab total energy in index",
            "ok" if db_energy is not None else "missing",
            f"{db_energy:.6f} eV" if db_energy is not None else "not recorded",
        ))
        # 4. Surface atoms indexed
        n_top = _count_rows(con, "top_atoms", surface_id)
        checks.append(_combine_check(
            "Surface atoms indexed",
            db_present=n_top > 0,
            disk_present=(
                (paths["metadata"] / "top_atoms.jsonl").is_file()
                and (paths["metadata"] / "atom_grid.npy").is_file()
            ),
            detail=f"{n_top} top atoms",
        ))
        # 5. Adsorption sites detected
        n_sites = _count_rows(con, "fcc_sites", surface_id)
        checks.append(_combine_check(
            "Adsorption sites detected",
            db_present=n_sites > 0,
            disk_present=(
                (paths["metadata"] / "fcc_sites.jsonl").is_file()
                and (paths["metadata"] / "site_grid.npy").is_file()
            ),
            detail=f"{n_sites} FCC sites",
        ))
        # 6. Hamiltonian matrix
        n_ham = _count_rows(con, "hamiltonian_exports", surface_id)
        checks.append(_combine_check(
            "Hamiltonian matrix stored",
            db_present=n_ham > 0,
            disk_present=(
                paths["openmx_slab"] / "hamiltonian_d_surface.npz"
            ).is_file(),
            detail=f"{n_ham} export(s)",
        ))

        core_complete = all(status_ == "ok" for _, status_, _ in checks)

        # Optional adsorbate coverage
        adsorbate_reports = []
        if args.adsorbates:
            site_ids = [
                r[0] for r in con.execute(
                    "SELECT site_id FROM fcc_sites WHERE surface_id=? "
                    "ORDER BY site_index", (surface_id,)
                )
            ] if con is not None else []
            for adsorbate in args.adsorbates:
                adsorbate_reports.append(_check_adsorbate_coverage(
                    con, paths, surface_id, adsorbate, site_ids
                ))
        adsorbates_complete = all(r["complete"] for r in adsorbate_reports)
        overall_complete = core_complete and (
            adsorbates_complete if adsorbate_reports else True
        )
    finally:
        if con is not None:
            con.close()

    if args.as_json:
        report = {
            "surface_id": surface_id,
            "status": status,
            "core_checks": [
                {"item": label, "status": st, "detail": detail}
                for label, st, detail in checks
            ],
            "core_complete": core_complete,
            "adsorbates": adsorbate_reports,
            "adsorbates_complete": adsorbates_complete if adsorbate_reports else None,
            "complete": overall_complete,
        }
        print(json.dumps(report, indent=2, sort_keys=True))
        return 0 if overall_complete else 1

    print("=" * 60)
    print(f"Surface : {surface_id}")
    print(f"Status  : {status or 'unknown'}")
    print("=" * 60)
    for label, st, detail in checks:
        line = f"{_CHECK_TAG[st]} {label}"
        if detail:
            line += f"  ({detail})"
        print(line)

    for rep in adsorbate_reports:
        print("-" * 60)
        total = rep["total_sites"]
        tag = "[COMPLETE]" if rep["complete"] else "[INCOMPLETE]"
        print(
            f"Adsorbate {rep['adsorbate']:<4s}: "
            f"configs {rep['registered']}/{total}  "
            f"structures {rep['structures']}/{total}  "
            f"energies {rep['energies']}/{total}   {tag}"
        )
        if total == 0:
            print("           no FCC sites to cover (run detect-sites first)")
            continue
        if rep["missing_config"]:
            print("           missing config   : "
                  + _fmt_id_list(rep["missing_config"], max_list))
        if rep["missing_structure"]:
            print("           missing structure: "
                  + _fmt_id_list(rep["missing_structure"], max_list))
        if rep["missing_energy"]:
            print("           missing energy   : "
                  + _fmt_id_list(rep["missing_energy"], max_list))

    print("=" * 60)
    if overall_complete:
        print("Result  : COMPLETE")
    else:
        parts = []
        if not core_complete:
            n_bad = sum(1 for _, st, _ in checks if st != "ok")
            parts.append(f"{n_bad} core item(s) missing/inconsistent")
        if adsorbate_reports and not adsorbates_complete:
            n_bad = sum(1 for r in adsorbate_reports if not r["complete"])
            parts.append(f"{n_bad} adsorbate(s) incomplete")
        print(f"Result  : INCOMPLETE ({'; '.join(parts)})")
    return 0 if overall_complete else 1


# ---------------------------------------------------------------------------
# Dataset merge
# ---------------------------------------------------------------------------

# Tables carrying per-surface rows, listed so that a row is always inserted
# after the row it references (surfaces first; adsorbate_configs after
# fcc_sites, which it points at through a composite foreign key).
_MERGE_TABLES = [
    "surfaces",
    "artifacts",
    "top_atoms",
    "fcc_sites",
    "adsorbate_configs",
    "hamiltonian_exports",
]

# Columns holding filesystem paths that live *inside* the surface directory and
# therefore must be re-rooted under the destination. Relative columns
# (surfaces.initial_cif/relaxed_cif) and external columns
# (hamiltonian_exports.scfout_path, which points outside the dataset) are left
# untouched on purpose.
_MERGE_PATH_COLS = {
    "surfaces": ["path"],
    "artifacts": ["path"],
    "adsorbate_configs": ["path"],
    "hamiltonian_exports": ["output_npz"],
}

# AUTOINCREMENT primary keys must be dropped on insert so the destination
# assigns fresh ids and two datasets never collide on them.
_MERGE_AUTOINC_COLS = {
    "artifacts": ["artifact_id"],
    "hamiltonian_exports": ["export_id"],
}


def _read_schema_version(db_path):
    with sqlite3.connect(db_path) as con:
        has_info = con.execute(
            "SELECT 1 FROM sqlite_master WHERE type='table' AND name='dataset_info'"
        ).fetchone()
        if not has_info:
            return None
        row = con.execute(
            "SELECT value FROM dataset_info WHERE key='schema_version'"
        ).fetchone()
        return int(row[0]) if row else None


def _rebase_surface_path(path_str, surface_id, dest_root):
    """Re-root a stored path under ``dest_root/surface_id``.

    Stored paths were written with whatever root the source dataset used (which
    may be relative or absolute, on another machine, etc.), but they always run
    through the ``<root>/<surface_id>/...`` layout. We locate the surface_id
    segment and rebuild everything after it beneath the destination, so the
    result is correct regardless of the original base. Paths without the
    surface_id segment (e.g. an external .scfout source) are returned unchanged.
    """
    if path_str is None:
        return None
    parts = Path(path_str).parts
    idx = None
    for i in range(len(parts) - 1, -1, -1):
        if parts[i] == surface_id:
            idx = i
            break
    if idx is None:
        return path_str
    new_path = Path(dest_root) / surface_id
    for segment in parts[idx + 1:]:
        new_path = new_path / segment
    return str(new_path)


def _prompt_conflict(surface_id):
    """Ask the user how to resolve a surface_id already present in the target."""
    print(f"[conflict] surface_id already exists in destination: {surface_id}")
    while True:
        answer = input("  [o]verwrite / [s]kip / [a]bort merge? ").strip().lower()
        if answer in ("o", "overwrite"):
            return "overwrite"
        if answer in ("s", "skip"):
            return "skip"
        if answer in ("a", "abort"):
            raise SystemExit("[abort] Merge cancelled by user.")
        print("  Please answer o, s or a.")


def _delete_surface_everywhere(con, dest_root, surface_id):
    """Remove a surface from the destination index rows and disk directory."""
    for table in _MERGE_TABLES:
        con.execute(f"DELETE FROM {table} WHERE surface_id=?", (surface_id,))
    dst_dir = Path(dest_root) / surface_id
    if dst_dir.is_dir():
        shutil.rmtree(dst_dir)


def _copy_surface_dir(src_root, dest_root, surface_id):
    src_dir = Path(src_root) / surface_id
    dst_dir = Path(dest_root) / surface_id
    if not src_dir.is_dir():
        print(f"[warn] No structure directory on disk for {surface_id}; "
              "copying index rows only.")
        return
    if dst_dir.exists():
        shutil.rmtree(dst_dir)
    shutil.copytree(src_dir, dst_dir)


def _copy_surface_rows(src_con, dst_con, surface_id, dest_root):
    """Copy every index row belonging to one surface, rewriting path columns."""
    for table in _MERGE_TABLES:
        cur = src_con.execute(
            f"SELECT * FROM {table} WHERE surface_id=?", (surface_id,)
        )
        columns = [desc[0] for desc in cur.description]
        for row in cur.fetchall():
            record = dict(zip(columns, row))
            for path_col in _MERGE_PATH_COLS.get(table, ()):
                if path_col in record:
                    record[path_col] = _rebase_surface_path(
                        record[path_col], surface_id, dest_root
                    )
            for drop_col in _MERGE_AUTOINC_COLS.get(table, ()):
                record.pop(drop_col, None)
            insert_cols = list(record.keys())
            placeholders = ", ".join("?" * len(insert_cols))
            dst_con.execute(
                f"INSERT INTO {table} ({', '.join(insert_cols)}) "
                f"VALUES ({placeholders})",
                [record[col] for col in insert_cols],
            )


def merge_datasets(args):
    """Merge every surface from --source into the --root dataset.

    Both datasets must share the current schema version. Each source surface is
    copied as its on-disk directory plus its index rows, with in-dataset path
    columns re-rooted under the destination. When a surface_id already exists in
    the destination the user is asked (by default) to overwrite, skip or abort.
    """
    src_root = Path(args.source)
    dst_root = Path(args.root)
    src_db = dataset_paths(src_root)["db"]

    if not src_db.is_file():
        raise SystemExit(f"[error] Source dataset index not found: {src_db}")
    if src_root.resolve() == dst_root.resolve():
        raise SystemExit("[error] Source and destination are the same dataset.")

    src_version = _read_schema_version(src_db)
    if src_version is None:
        raise SystemExit(
            f"[error] Source dataset has no schema_version: {src_root}"
        )
    if src_version != SCHEMA_VERSION:
        raise SystemExit(
            f"[error] Source schema v{src_version} is incompatible with "
            f"v{SCHEMA_VERSION}."
        )

    # Creates the destination if absent and validates its schema otherwise.
    init_db(dst_root)

    src_con = sqlite3.connect(src_db)
    dst_con = sqlite3.connect(dataset_paths(dst_root)["db"])
    try:
        src_ids = [
            row[0] for row in src_con.execute(
                "SELECT surface_id FROM surfaces ORDER BY surface_id"
            )
        ]
        if not src_ids:
            print("[info] Source dataset has no surfaces; nothing to merge.")
            return

        added = overwritten = skipped = 0
        for surface_id in src_ids:
            if surface_exists(dst_root, surface_id):
                action = args.on_conflict
                if action == "ask":
                    action = _prompt_conflict(surface_id)
                if action == "skip":
                    print(f"[skip] {surface_id}")
                    skipped += 1
                    continue
                _delete_surface_everywhere(dst_con, dst_root, surface_id)
                overwritten += 1
            else:
                added += 1

            _copy_surface_dir(src_root, dst_root, surface_id)
            _copy_surface_rows(src_con, dst_con, surface_id, dst_root)
            dst_con.commit()
            print(f"[merge] {surface_id}")

        print("=" * 60)
        print(f"[done] Merge complete: {added} added, {overwritten} overwritten, "
              f"{skipped} skipped.")
    finally:
        src_con.close()
        dst_con.close()


def delete_surface(args):
    """Delete one surface: its index rows across all tables and its directory.

    Rows and the on-disk directory are two independent stores, so both are
    removed together to avoid leaving orphan rows or an orphan folder behind.
    Deletion is irreversible (shutil.rmtree), hence the confirmation prompt
    unless --yes is given.
    """
    surface_id = args.surface_id
    if not surface_exists(args.root, surface_id):
        raise SystemExit(f"[error] Surface not found in dataset: {surface_id}")

    paths = dataset_paths(args.root, surface_id)
    surface_dir = paths["surface"]

    if not args.yes:
        print(f"About to permanently delete surface: {surface_id}")
        print(f"  index rows in : {', '.join(_MERGE_TABLES)}")
        print(f"  directory     : {surface_dir}")
        answer = input("Proceed? [y/N] ").strip().lower()
        if answer not in ("y", "yes"):
            raise SystemExit("[abort] Deletion cancelled.")

    db_path = dataset_paths(args.root)["db"]
    if db_path.is_file():
        # _delete_surface_everywhere removes both the index rows and the
        # on-disk directory for this surface_id.
        with sqlite3.connect(db_path) as con:
            _delete_surface_everywhere(con, args.root, surface_id)
            con.commit()
    elif surface_dir.is_dir():
        # No index present, but a structure directory exists on disk.
        shutil.rmtree(surface_dir)

    print(f"[done] Deleted surface {surface_id}")


def add_common_surface_args(p):
    p.add_argument("--root", default="dataset", help="Dataset root directory.")
    p.add_argument("--surface-id", required=True, help="Composition-based surface ID.")


def build_parser():
    p = argparse.ArgumentParser(
        description="Manage HEA surface dataset layout, metadata and SQLite index.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    sub = p.add_subparsers(dest="cmd", required=True)

    sp = sub.add_parser("init", help="Initialize dataset directory and SQLite index.")
    sp.add_argument("--root", default="dataset")
    sp.set_defaults(func=lambda args: init_db(args.root))

    sp = sub.add_parser(
        "create-sample",
        help="Generate an SQS slab from elements (+optional ratios) and register it.",
    )
    sp.add_argument("--root", default="dataset")
    sp.add_argument(
        "-e", "--elements", nargs="+", required=True, metavar="EL",
        help="Element types to include, e.g. Fe Co Ni Cr Mn.",
    )
    sp.add_argument(
        "-r", "--ratios", nargs="+", type=float, default=None, metavar="R",
        help="Ratio per element. Randomly generated (and retried on collision) "
             "when omitted.",
    )
    sp.add_argument(
        "-s", "--size", nargs=3, type=int, default=(4, 4, 4),
        metavar=("NX", "NY", "NZ"), help="Repeats/layers along x, y and z.",
    )
    sp.add_argument(
        "--vacuum", type=float, default=15.0,
        help="Total vacuum-layer thickness (Angstrom).",
    )
    sp.add_argument(
        "-a", "--lattice-constant", type=float, default=None,
        help="Template lattice constant; defaults to a Vegard's-law estimate.",
    )
    sp.add_argument(
        "--cutoffs", nargs="+", type=float, default=[6.0, 4.5],
        help="icet cluster-space cutoff radii (Angstrom): pair, triplet, ...",
    )
    sp.add_argument(
        "--n-steps", type=int, default=10000,
        help="Monte Carlo steps for SQS optimization.",
    )
    sp.add_argument("--seed", type=int, default=None, help="Random seed.")
    sp.add_argument(
        "--no-sqs", action="store_true",
        help="Skip SQS and register the randomly substituted structure.",
    )
    sp.add_argument(
        "--max-attempts", type=int, default=DEFAULT_MAX_ATTEMPTS,
        help="Max random-ratio draws before giving up on a novel composition.",
    )
    sp.set_defaults(func=create_sample)

    sp = sub.add_parser(
        "record-relaxed",
        help="Register an externally relaxed slab for an existing composition.",
    )
    add_common_surface_args(sp)
    sp.add_argument("--relaxed-cif", required=True, help="Relaxed slab CIF.")
    sp.add_argument(
        "--energy", type=float, default=None,
        help="Total DFT energy of the relaxed slab in eV (from OUTCAR TOTEN).",
    )
    sp.set_defaults(func=record_relaxed)

    sp = sub.add_parser("index-surface", help="Build top atom metadata/grid.")
    add_common_surface_args(sp)
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
    add_common_surface_args(sp)
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
        help="Create per-site records and initial adsorbate CIF structures.",
    )
    add_common_surface_args(sp)
    sp.add_argument(
        "--adsorbate", required=True, choices=["N", "NH", "NH2", "NH3"],
        help="Adsorbate placed with N bound to the surface.",
    )
    sp.add_argument("--side", choices=["top", "bottom"], default="top")
    sp.add_argument("--height", type=float, default=1.25)
    sp.add_argument("--nh", type=float, default=1.02, help="N-H bond length.")
    sp.add_argument(
        "--rebuild-initial", action="store_true",
        help="Regenerate the initial adsorbate structure and reset the energy "
             "status even for sites that already have a record. Off by default "
             "so reruns never discard computed adsorption energies.",
    )
    sp.set_defaults(func=create_adsorbate_records)

    sp = sub.add_parser("record-energy", help="Record a manual adsorption energy.")
    add_common_surface_args(sp)
    sp.add_argument("--adsorbate", required=True)
    sp.add_argument(
        "--site-id", type=int, required=True,
        help="Positive site number, e.g. 1 for site_0001.",
    )
    sp.add_argument("--relaxed-cif", required=True, help="Relaxed adsorbate CIF.")
    sp.add_argument("--energy", type=float, required=True)
    sp.add_argument(
        "--max-adsorbate-displacement", type=float, default=2.0,
        help="Maximum allowed adsorbate-atom displacement in Angstrom.",
    )
    sp.add_argument("--status", default="manually_entered")
    sp.add_argument("--notes", default="")
    sp.set_defaults(func=record_energy)

    sp = sub.add_parser(
        "extract-hamiltonian",
        help="Extract top-layer d-orbital Hamiltonian data from OpenMX output.",
    )
    add_common_surface_args(sp)
    sp.add_argument(
        "--scfout", required=True, help="OpenMX .scfout file.",
    )
    sp.add_argument(
        "--dat", required=True, help="OpenMX .dat input used for the run.",
    )
    sp.add_argument("-o", "--output", default=None)
    sp.add_argument("--basis-output", default=None)
    sp.set_defaults(func=extract_hamiltonian)

    sp = sub.add_parser(
        "check",
        help="Check completeness of a surface (structures, energy, indexing, "
             "Hamiltonian) and optionally per-adsorbate site coverage.",
    )
    add_common_surface_args(sp)
    sp.add_argument(
        "-a", "--adsorbates", nargs="+", default=None, metavar="SPECIES",
        help="Also verify per-site coverage for these adsorbate species, e.g. "
             "N NH NH2. Each site must have a registered config, a relaxed "
             "structure on disk and an adsorption energy in the index.",
    )
    sp.add_argument(
        "--max-list", type=int, default=10,
        help="Max site IDs listed per missing category (0 = list all).",
    )
    sp.add_argument(
        "--json", dest="as_json", action="store_true",
        help="Emit the report as JSON instead of a human-readable table.",
    )
    sp.set_defaults(func=check_surface)

    sp = sub.add_parser(
        "delete",
        help="Delete a surface: its index rows and its on-disk directory.",
    )
    add_common_surface_args(sp)
    sp.add_argument(
        "-y", "--yes", action="store_true",
        help="Skip the confirmation prompt (deletion is irreversible).",
    )
    sp.set_defaults(func=delete_surface)

    sp = sub.add_parser(
        "merge",
        help="Merge all surfaces from a source dataset into --root.",
    )
    sp.add_argument("--root", default="dataset", help="Destination dataset root.")
    sp.add_argument(
        "--source", required=True, help="Source dataset root to merge from.",
    )
    sp.add_argument(
        "--on-conflict", choices=["ask", "skip", "overwrite"], default="ask",
        help="How to resolve a surface_id that already exists in the "
             "destination. 'ask' prompts on the command line for each conflict.",
    )
    sp.set_defaults(func=merge_datasets)

    return p


def main(argv=None):
    args = build_parser().parse_args(argv)
    rc = args.func(args)
    # Commands that report a status (e.g. check) return an int exit code;
    # the rest return None and are treated as success.
    return rc if isinstance(rc, int) else 0


if __name__ == "__main__":
    raise SystemExit(main())
