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

## 2. Register the Clean Slab

```bash
python hea_dataset.py create-sample \
  --root dataset \
  --initial-cif test_sqs.cif \
  --relaxed-cif test_sqs-opt.cif
```

Both CIF files must have identical element counts. For the example composition,
the resulting paths are:

```text
dataset/Co_13-Cr_13-Fe_13-Mn_12-Ni_13/manifest.json
dataset/Co_13-Cr_13-Fe_13-Mn_12-Ni_13/structures/00_initial_sqs.cif
dataset/Co_13-Cr_13-Fe_13-Mn_12-Ni_13/structures/01_relaxed_slab.cif
dataset/Co_13-Cr_13-Fe_13-Mn_12-Ni_13/metadata/
dataset/Co_13-Cr_13-Fe_13-Mn_12-Ni_13/adsorbates/
dataset/Co_13-Cr_13-Fe_13-Mn_12-Ni_13/openmx_slab/
```

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
