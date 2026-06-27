# HEA Dataset Workflow

This repository stores one composition as one `surface_id`:

```text
Fe_20-Co_20-Ni_20-Cr_20-Mn_20
```

Different random/SQS structures under the same composition are separated by
`sample_id`, for example `sample_0001`.

## 1. Initialize Dataset

```bash
python hea_dataset.py init --root dataset
```

## 2. Create A Sample

```bash
python hea_dataset.py create-sample \
  --root dataset \
  --composition Fe=20 Co=20 Ni=20 Cr=20 Mn=20 \
  --sample-id sample_0001 \
  --initial-cif HEA_FeCoNiCrMn_111_sqs.cif \
  --relaxed-cif HEA_FeCoNiCrMn_111_relaxed.cif
```

This creates:

```text
dataset/surfaces/Fe_20-Co_20-Ni_20-Cr_20-Mn_20/sample_0001/
  manifest.json
  structures/
    00_initial_sqs.cif
    01_relaxed_slab.cif
  metadata/
  openmx_slab/
  adsorbates/
```

## 3. Index Surface Atoms

```bash
python hea_dataset.py index-surface \
  --root dataset \
  --surface-id Fe_20-Co_20-Ni_20-Cr_20-Mn_20 \
  --sample-id sample_0001
```

Outputs:

```text
metadata/top_atoms.jsonl
metadata/atom_grid.npy
```

The grid is sorted by wrapped fractional `y`, then fractional `x`.

## 4. Detect FCC Sites

```bash
python hea_dataset.py detect-sites \
  --root dataset \
  --surface-id Fe_20-Co_20-Ni_20-Cr_20-Mn_20 \
  --sample-id sample_0001
```

Outputs:

```text
metadata/fcc_sites.jsonl
metadata/site_grid.npy
```

## 5. Create Adsorbate Records

```bash
python hea_dataset.py create-adsorbate-records \
  --root dataset \
  --surface-id Fe_20-Co_20-Ni_20-Cr_20-Mn_20 \
  --sample-id sample_0001 \
  --adsorbate N
```

Each site folder contains only:

```text
adsorbates/N/site_0001/
  00_initial_adsorbate.cif
  01_relaxed_adsorbate.cif
  adsorption_energy.json
```

The two CIF files are placeholders for the structures generated/relaxed by your
external workflow. The JSON file is used for manual adsorption energy input.

## 6. Record Adsorption Energy

```bash
python hea_dataset.py record-energy \
  --root dataset \
  --surface-id Fe_20-Co_20-Ni_20-Cr_20-Mn_20 \
  --sample-id sample_0001 \
  --adsorbate N \
  --site-id site_0001 \
  --energy -1.23
```

## 7. Extract Surface d-Orbital Hamiltonian

After `cif_to_openmx.py` generates the clean slab OpenMX input and OpenMX writes
`result.scfout`, run:

```bash
python extract_openmx_hamiltonian.py \
  --scfout dataset/surfaces/Fe_20-Co_20-Ni_20-Cr_20-Mn_20/sample_0001/openmx_slab/result.scfout \
  --dat dataset/surfaces/Fe_20-Co_20-Ni_20-Cr_20-Mn_20/sample_0001/openmx_slab/input.dat \
  --top-atoms dataset/surfaces/Fe_20-Co_20-Ni_20-Cr_20-Mn_20/sample_0001/metadata/top_atoms.jsonl \
  --output dataset/surfaces/Fe_20-Co_20-Ni_20-Cr_20-Mn_20/sample_0001/openmx_slab/hamiltonian_d_surface.npz \
  --dataset-root dataset \
  --surface-id Fe_20-Co_20-Ni_20-Cr_20-Mn_20 \
  --sample-id sample_0001
```

The `.npz` contains:

```text
H_d
d_basis_indices
d_labels
surface_atom_ids
openmx_atom_indices
spin_switch
source_scfout
source_dat
```
