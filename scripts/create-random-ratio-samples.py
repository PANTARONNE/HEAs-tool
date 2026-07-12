#!/usr/bin/env python3
"""Create a parallel batch of screened, random-ratio quinary HEA samples.

Compositions are selected centrally so they are unique, then generated in
isolated dataset shards.  Completed shards are merged serially into the target
dataset, avoiding concurrent writes to the target SQLite database.
"""

from __future__ import annotations

import argparse
from concurrent.futures import ProcessPoolExecutor, as_completed
from contextlib import redirect_stderr, redirect_stdout
import json
import os
from pathlib import Path
import shutil
import sys
from types import SimpleNamespace
import uuid

import numpy as np


PROJECT_ROOT = Path(__file__).resolve().parents[1]
if str(PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(PROJECT_ROOT))

import build_hea_surface as builder  # noqa: E402
import hea_dataset  # noqa: E402


ELEMENT_POOL = ("Fe", "Co", "Ni", "Cu", "Mo", "Zn", "Ga", "In", "Sn", "W")
CRITERION_KEYS = ("entropy_ok", "size_ok", "hmix_ok", "vec_ok")


def parse_args(argv=None):
    parser = argparse.ArgumentParser(
        description=(
            "Create screened random-ratio quinary HEA structures in parallel "
            "by calling hea_dataset.create_sample."
        ),
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument("--root", default="randomRatio", help="Output dataset root.")
    parser.add_argument(
        "--count", type=int, default=2000,
        help="Target total number of matching samples in the dataset (resume-safe).",
    )
    parser.add_argument(
        "--workers", type=int, default=min(4, os.cpu_count() or 1),
        help="Parallel worker processes. SQS generation can use substantial memory.",
    )
    parser.add_argument(
        "--chunk-size", type=int, default=25,
        help="Samples per temporary dataset shard before it is merged.",
    )
    parser.add_argument("--seed", type=int, default=20260713, help="Batch random seed.")
    parser.add_argument(
        "-s", "--size", nargs=3, type=int, default=(4, 4, 4),
        metavar=("NX", "NY", "NZ"), help="FCC(111) slab repeats/layers.",
    )
    parser.add_argument("--vacuum", type=float, default=15.0)
    parser.add_argument("--cutoffs", nargs="+", type=float, default=[6.0, 4.5])
    parser.add_argument("--n-steps", type=int, default=10000)
    parser.add_argument(
        "--no-sqs", action="store_true",
        help="Skip SQS optimization and retain random site substitution.",
    )
    parser.add_argument(
        "--max-draws", type=int, default=2_000_000,
        help="Maximum candidate compositions examined during screening.",
    )
    return parser.parse_args(argv)


def matching_existing_ids(root: Path) -> set[str]:
    """Return registered pool-derived quinary samples passing all criteria."""
    db = hea_dataset.dataset_paths(root)["db"]
    if not db.is_file():
        return set()

    import sqlite3

    result = set()
    with sqlite3.connect(db) as connection:
        table = connection.execute(
            "SELECT 1 FROM sqlite_master WHERE type='table' AND name='surfaces'"
        ).fetchone()
        if not table:
            return set()
        for surface_id, composition_json in connection.execute(
            "SELECT surface_id, composition_json FROM surfaces"
        ):
            composition = json.loads(composition_json)
            elements = [row["element"] for row in composition]
            counts = np.asarray(
                [row["atom_count"] for row in composition], dtype=float
            )
            if len(set(elements)) != 5 or not set(elements).issubset(ELEMENT_POOL):
                continue
            passed, report = hea_dataset._check_hea_criteria(
                elements, counts / counts.sum()
            )
            if passed and all(report[key] is True for key in CRITERION_KEYS):
                result.add(surface_id)
    return result


def select_compositions(args, excluded_ids: set[str], needed: int):
    """Draw unique integer compositions that pass all four HEA checks."""
    rng = np.random.default_rng(args.seed)
    template = builder.build_template_slab(
        ELEMENT_POOL[0], tuple(args.size), builder.DEFAULT_A, args.vacuum
    )
    n_sites = len(template)
    selected = []
    seen = set(excluded_ids)

    for draw in range(1, args.max_draws + 1):
        # Randomly choose exactly five distinct elements and random positive ratios.
        elements = list(rng.choice(ELEMENT_POOL, size=5, replace=False))
        fractions = rng.random(5) + 0.1
        fractions /= fractions.sum()
        counts = builder.largest_remainder_counts(fractions, n_sites)
        if np.any(counts <= 0):
            continue

        # create-sample screens its supplied ratios before constructing the slab.
        # Screen the exact integer ratios here so its result cannot differ due to
        # largest-remainder rounding.
        exact_fractions = counts / counts.sum()
        passed, report = hea_dataset._check_hea_criteria(elements, exact_fractions)
        if not passed or not all(report[key] is True for key in CRITERION_KEYS):
            continue

        composition = hea_dataset.composition_from_counts(elements, counts)
        surface_id = hea_dataset.surface_id_from_composition(composition)
        if surface_id in seen:
            continue
        seen.add(surface_id)
        selected.append({
            "surface_id": surface_id,
            "elements": elements,
            "counts": [int(value) for value in counts],
            "seed": int(rng.integers(0, 2**31 - 1)),
            "criteria": report,
        })
        if len(selected) == needed:
            return selected, n_sites, draw

    raise RuntimeError(
        f"Only selected {len(selected)} of {needed} required compositions after "
        f"{args.max_draws} draws. Increase --max-draws or change the slab size."
    )


def generate_shard(payload):
    """Generate one isolated shard; designed to run in a child process."""
    shard_root = Path(payload["shard_root"])
    samples = payload["samples"]
    common = payload["common"]
    shard_root.mkdir(parents=True, exist_ok=False)
    log_path = shard_root / "generation.log"

    # Avoid multiplying native-library threads inside each worker process.
    os.environ.setdefault("OMP_NUM_THREADS", "1")
    with log_path.open("a", encoding="utf-8") as log, \
            redirect_stdout(log), redirect_stderr(log):
        for sample in samples:
            call_args = SimpleNamespace(
                root=str(shard_root),
                elements=sample["elements"],
                ratios=sample["counts"],
                size=common["size"],
                vacuum=common["vacuum"],
                lattice_constant=None,
                cutoffs=common["cutoffs"],
                n_steps=common["n_steps"],
                seed=sample["seed"],
                no_sqs=common["no_sqs"],
                max_attempts=1,
            )
            hea_dataset.create_sample(call_args)
    return str(shard_root), len(samples)


def merge_shard(shard_root: Path, target_root: Path):
    args = SimpleNamespace(
        source=str(shard_root), root=str(target_root), on_conflict="skip"
    )
    # Batch progress is clearer than thousands of per-surface merge lines.
    with (shard_root / "merge.log").open("w", encoding="utf-8") as log, \
            redirect_stdout(log), redirect_stderr(log):
        hea_dataset.merge_datasets(args)


def write_plan(path: Path, args, selected, n_sites, draws, existing_count):
    path.parent.mkdir(parents=True, exist_ok=True)
    data = {
        "element_pool": list(ELEMENT_POOL),
        "target_count": args.count,
        "existing_count": existing_count,
        "new_count": len(selected),
        "seed": args.seed,
        "size": list(args.size),
        "n_sites": n_sites,
        "candidate_draws": draws,
        "samples": selected,
    }
    path.write_text(json.dumps(data, indent=2, sort_keys=True) + "\n", encoding="utf-8")


def main(argv=None):
    args = parse_args(argv)
    if args.count < 1 or args.workers < 1 or args.chunk_size < 1:
        raise SystemExit("[error] --count, --workers and --chunk-size must be positive.")

    target_root = Path(args.root).resolve()
    existing_ids = matching_existing_ids(target_root)
    needed = max(0, args.count - len(existing_ids))
    if needed == 0:
        print(f"[done] {target_root} already contains {len(existing_ids)} matching samples.")
        return 0

    print(f"[select] Need {needed} new samples ({len(existing_ids)} already present).")
    selected, n_sites, draws = select_compositions(args, existing_ids, needed)
    print(f"[select] Accepted {needed} unique compositions from {draws} draws.")

    run_id = uuid.uuid4().hex[:10]
    shard_parent = target_root / f".batch_shards_{run_id}"
    shard_parent.mkdir(parents=True, exist_ok=False)
    write_plan(
        shard_parent / "batch_plan.json", args, selected, n_sites, draws,
        len(existing_ids),
    )

    chunks = [selected[i:i + args.chunk_size]
              for i in range(0, len(selected), args.chunk_size)]
    common = {
        "size": tuple(args.size), "vacuum": args.vacuum,
        "cutoffs": list(args.cutoffs), "n_steps": args.n_steps,
        "no_sqs": args.no_sqs,
    }
    payloads = [{
        "shard_root": str(shard_parent / f"shard_{index:05d}"),
        "samples": chunk,
        "common": common,
    } for index, chunk in enumerate(chunks, start=1)]

    completed = 0
    try:
        with ProcessPoolExecutor(max_workers=args.workers) as executor:
            futures = {executor.submit(generate_shard, payload): payload
                       for payload in payloads}
            for future in as_completed(futures):
                shard_root_text, shard_count = future.result()
                shard_root = Path(shard_root_text)
                merge_shard(shard_root, target_root)
                completed += shard_count
                print(f"[progress] {completed}/{needed} new samples merged")
                shutil.rmtree(shard_root)
    except BaseException:
        print(
            f"[error] Batch stopped after merging {completed}/{needed} samples. "
            f"Logs and unmerged shards remain in {shard_parent}", file=sys.stderr,
        )
        raise

    # Keep the reproducibility record in the completed dataset.
    plan_destination = target_root / f"batch_plan_{run_id}.json"
    shutil.move(str(shard_parent / "batch_plan.json"), plan_destination)
    shutil.rmtree(shard_parent)
    print(f"[done] Dataset now contains {len(existing_ids) + completed} matching samples.")
    print(f"[done] Batch plan: {plan_destination}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
