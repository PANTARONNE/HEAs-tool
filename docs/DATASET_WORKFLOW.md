# HEA Dataset Workflow

The dataset uses one directory per composition. The directory name is derived
from the CIF atom counts with element symbols sorted alphabetically. For
example, a structure containing Co13Cr13Fe13Mn12Ni13 is stored as:

```text
dataset/Co_13-Cr_13-Fe_13-Mn_12-Ni_13/
```

Only one structure is allowed for each composition. `create-sample` rejects a
composition already present in SQLite or on disk.

## 1. Initialize the Dataset

```bash
python hea_dataset.py init --root dataset
```

This creates the dataset-level files:

```text
dataset/index.sqlite
dataset/dataset_manifest.json
```

Calling `create-sample` also initializes these files when necessary.

## 2. Generate and Register the Clean Slab

`create-sample` builds the FCC(111) SQS slab directly from elements and optional
ratios. When `--ratios` is omitted, ratios are drawn at random; if the resulting
composition already exists, new ratios are drawn until a novel one is found
(bounded by `--max-attempts`). Fixed ratios that collide are a hard error.

```bash
python hea_dataset.py create-sample \
  --root dataset \
  --elements Fe Co Ni Cr Mn \
  --ratios 1 1 1 1 1
```

Only the initial SQS structure is registered here; the relaxed slab is added
later with `record-relaxed`. The resulting paths are:

```text
dataset/Co_13-Cr_13-Fe_13-Mn_12-Ni_13/manifest.json
dataset/Co_13-Cr_13-Fe_13-Mn_12-Ni_13/structures/00_initial_sqs.cif
dataset/Co_13-Cr_13-Fe_13-Mn_12-Ni_13/metadata/
dataset/Co_13-Cr_13-Fe_13-Mn_12-Ni_13/adsorbates/
dataset/Co_13-Cr_13-Fe_13-Mn_12-Ni_13/openmx_slab/
```

Omit `--ratios` to let a random composition be chosen. Use `--no-sqs` to skip
SQS optimization, and `--seed` for reproducibility.

## 2b. Register the Relaxed Slab

After relaxing `00_initial_sqs.cif` externally, register the result. The
command checks the composition matches and that atom count and order are
preserved (index-surface associates relaxed coordinates by ASE atom index):

```bash
python hea_dataset.py record-relaxed \
  --root dataset \
  --surface-id Co_13-Cr_13-Fe_13-Mn_12-Ni_13 \
  --relaxed-cif test_sqs-opt.cif
```

This copies the structure to `structures/01_relaxed_slab.cif` and updates the
manifest and SQLite index (status `slab_relaxed`).

## 3. Index Surface Atoms

```bash
python hea_dataset.py index-surface \
  --root dataset \
  --surface-id Co_13-Cr_13-Fe_13-Mn_12-Ni_13
```

By default this reads the two registered CIF files and writes:

```text
dataset/Co_13-Cr_13-Fe_13-Mn_12-Ni_13/metadata/top_atoms.jsonl
dataset/Co_13-Cr_13-Fe_13-Mn_12-Ni_13/metadata/atom_grid.npy
```

The initial CIF defines atom IDs and the surface grid. The relaxed coordinates
are associated by ASE atom index, so both CIF files must preserve atom order.

## 4. Detect FCC Sites

```bash
python hea_dataset.py detect-sites \
  --root dataset \
  --surface-id Co_13-Cr_13-Fe_13-Mn_12-Ni_13
```

The default input is the registered relaxed slab
`structures/01_relaxed_slab.cif`. Outputs:

```text
dataset/Co_13-Cr_13-Fe_13-Mn_12-Ni_13/metadata/fcc_sites.jsonl
dataset/Co_13-Cr_13-Fe_13-Mn_12-Ni_13/metadata/site_grid.npy
```

## 5. Create Adsorbate Structures and Records

```bash
python hea_dataset.py create-adsorbate-records \
  --root dataset \
  --surface-id Co_13-Cr_13-Fe_13-Mn_12-Ni_13 \
  --adsorbate N
```

Supported adsorbates are `N`, `NH`, `NH2`, and `NH3`.
For every recorded FCC site, the command reads the registered relaxed slab and creates:

```text
dataset/Co_13-Cr_13-Fe_13-Mn_12-Ni_13/adsorbates/N/site_0001/00_initial_adsorbate.cif
dataset/Co_13-Cr_13-Fe_13-Mn_12-Ni_13/adsorbates/N/site_0001/adsorption_energy.json
...
dataset/Co_13-Cr_13-Fe_13-Mn_12-Ni_13/adsorbates/N/site_0016/00_initial_adsorbate.cif
dataset/Co_13-Cr_13-Fe_13-Mn_12-Ni_13/adsorbates/N/site_0016/adsorption_energy.json
```

After external relaxation, place each result at:

```text
dataset/Co_13-Cr_13-Fe_13-Mn_12-Ni_13/adsorbates/N/site_0001/01_relaxed_adsorbate.cif
```

## 6. Record Adsorption Energy

```bash
python hea_dataset.py record-energy \
  --root dataset \
  --surface-id Co_13-Cr_13-Fe_13-Mn_12-Ni_13 \
  --adsorbate N \
  --site-id 1 \
  --relaxed-cif N-site-0001-relaxed.cif \
  --energy -1.23
```

Before recording, the command checks atom count and atom order, removes the
rigid translation measured from substrate atoms, and calculates minimum-image
adsorbate displacements. The default maximum allowed displacement is 2.0 Å;
override it with `--max-adsorbate-displacement`. A failed check does not copy
the CIF or update the energy. A successful check copies the structure to:

```text
dataset/Co_13-Cr_13-Fe_13-Mn_12-Ni_13/adsorbates/N/site_0001/01_relaxed_adsorbate.cif
```

It also writes the energy and displacement metrics to `adsorption_energy.json`
and updates the SQLite index.

## 7. Generate the Clean-Slab OpenMX Input

```bash
python cif_to_openmx.py \
  dataset/Co_13-Cr_13-Fe_13-Mn_12-Ni_13/structures/01_relaxed_slab.cif \
  --data-path DFT_DATA19 \
  -o dataset/Co_13-Cr_13-Fe_13-Mn_12-Ni_13/openmx_slab/input.dat
```

Run OpenMX externally and place its output at, for example:

```text
dataset/Co_13-Cr_13-Fe_13-Mn_12-Ni_13/openmx_slab/input.scfout
```

## 8. Extract the Surface d-Orbital Hamiltonian

```bash
python hea_dataset.py extract-hamiltonian \
  --root dataset \
  --surface-id Co_13-Cr_13-Fe_13-Mn_12-Ni_13 \
  --dat dataset/Co_13-Cr_13-Fe_13-Mn_12-Ni_13/openmx_slab/input.dat \
  --scfout dataset/Co_13-Cr_13-Fe_13-Mn_12-Ni_13/openmx_slab/input.scfout
```

The command obtains `metadata/top_atoms.jsonl` from `surface-id`
automatically. `--dat` and `--scfout` are required inputs. Use `--output` only
when a non-default output path is needed.

The basis mapping is written by default to:

```text
dataset/Co_13-Cr_13-Fe_13-Mn_12-Ni_13/openmx_slab/hamiltonian_d_surface.npz.basis.jsonl
```

## 9. Merge Two Datasets

`merge` combines every surface from a source dataset into a destination dataset.
Because a dataset is an SQLite index *plus* the per-surface directories it points
at, this cannot be done by copying rows alone: `merge` copies each
`<surface_id>/` directory to the destination and re-roots the stored path
columns so they point at their new location.

```bash
python hea_dataset.py merge \
  --root dataset \
  --source other_dataset
```

`--root` is the destination (rows are written into it); `--source` is the
dataset being merged in. Both datasets must share the same schema version, or
the command aborts before making any change.

For each source surface the command copies its directory and index rows, then
rewrites the in-dataset path columns (`surfaces.path`, `artifacts.path`,
`adsorbate_configs.path`, `hamiltonian_exports.output_npz`) to the destination
root. Relative columns (`initial_cif`, `relaxed_cif`) and external paths
(`hamiltonian_exports.scfout_path`) are left unchanged.

### Handling `surface_id` conflicts

Because `surface_id` is the composition-based primary key, a surface already
present in the destination is a conflict. By default `merge` prints the
conflict and prompts for a decision per surface:

```text
[conflict] surface_id already exists in destination: Co_13-Cr_13-Fe_13-Mn_12-Ni_13
  [o]verwrite / [s]kip / [a]bort merge?
```

- `o` — delete the destination surface (index rows and directory) and replace it
  with the source copy.
- `s` — leave the destination surface untouched and skip the source one.
- `a` — abort the whole merge.

For non-interactive runs, set the policy up front with `--on-conflict`:

```bash
python hea_dataset.py merge --root dataset --source other_dataset --on-conflict skip
```

Choices are `ask` (default), `skip`, and `overwrite`. The merge writes into the
destination in place, and `overwrite` deletes the replaced directory, so back up
the destination first if the source is not fully trusted.
