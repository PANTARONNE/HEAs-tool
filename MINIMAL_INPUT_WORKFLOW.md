# Minimal Input Workflow

This document lists the minimum manual inputs needed to create one complete HEA
training sample.

## Required Manual Inputs

```text
composition
sample_id
initial_slab_cif
relaxed_slab_cif
adsorbate
adsorption_energy_eV for each site
openmx_data_path
openmx_scfout
```

Example values:

```text
composition         Fe=20 Co=20 Ni=20 Cr=20 Mn=20
sample_id           sample_0001
initial_slab_cif    test_sqs.cif
relaxed_slab_cif    test_sqs-opt.cif
adsorbate           N
openmx_data_path    DFT_DATA19
openmx_scfout       openmx_slab/input.scfout
```

The `surface_id` is derived from `composition`:

```text
Fe_20-Co_20-Ni_20-Cr_20-Mn_20
```

## 1. Create Dataset Sample

```powershell
python hea_dataset.py create-sample --root dataset --composition Fe=20 Co=20 Ni=20 Cr=20 Mn=20 --sample-id sample_0001 --initial-cif test_sqs.cif --relaxed-cif test_sqs-opt.cif
```

This copies the CIF files to:

```text
dataset/surfaces/Fe_20-Co_20-Ni_20-Cr_20-Mn_20/sample_0001/structures/00_initial_sqs.cif
dataset/surfaces/Fe_20-Co_20-Ni_20-Cr_20-Mn_20/sample_0001/structures/01_relaxed_slab.cif
```

## 2. Index Surface Atoms

```powershell
python hea_dataset.py index-surface --root dataset --surface-id Fe_20-Co_20-Ni_20-Cr_20-Mn_20 --sample-id sample_0001
```

Outputs:

```text
metadata/top_atoms.jsonl
metadata/atom_grid.npy
```

## 3. Detect FCC Sites

```powershell
python hea_dataset.py detect-sites --root dataset --surface-id Fe_20-Co_20-Ni_20-Cr_20-Mn_20 --sample-id sample_0001
```

Outputs:

```text
metadata/fcc_sites.jsonl
metadata/site_grid.npy
```

## 4. Create Adsorbate Records

```powershell
python hea_dataset.py create-adsorbate-records --root dataset --surface-id Fe_20-Co_20-Ni_20-Cr_20-Mn_20 --sample-id sample_0001 --adsorbate N
```

This creates one folder per site:

```text
adsorbates/N/site_0001/
  adsorption_energy.json
```

The only expected files in each final adsorbate folder are:

```text
00_initial_adsorbate.cif
01_relaxed_adsorbate.cif
adsorption_energy.json
```

## 5. Generate Initial Adsorbate CIFs

```powershell
python add_fcc_adsorbate.py dataset/surfaces/Fe_20-Co_20-Ni_20-Cr_20-Mn_20/sample_0001/structures/01_relaxed_slab.cif --adsorbate N --site all
```

Move or rename the generated files into the corresponding site folders as:

```text
adsorbates/N/site_0001/00_initial_adsorbate.cif
adsorbates/N/site_0002/00_initial_adsorbate.cif
...
```

After external relaxation, place each relaxed adsorbate structure as:

```text
adsorbates/N/site_0001/01_relaxed_adsorbate.cif
adsorbates/N/site_0002/01_relaxed_adsorbate.cif
...
```

## 6. Record Adsorption Energies

```powershell
python hea_dataset.py record-energy --root dataset --surface-id Fe_20-Co_20-Ni_20-Cr_20-Mn_20 --sample-id sample_0001 --adsorbate N --site-id site_0001 --energy -1.23
```

Repeat for each site.

## 7. Generate OpenMX Input For Clean Slab

```powershell
python cif_to_openmx.py dataset/surfaces/Fe_20-Co_20-Ni_20-Cr_20-Mn_20/sample_0001/structures/01_relaxed_slab.cif --data-path DFT_DATA19 -o dataset/surfaces/Fe_20-Co_20-Ni_20-Cr_20-Mn_20/sample_0001/openmx_slab/input.dat
```

Run OpenMX externally and place the resulting `.scfout` in:

```text
dataset/surfaces/Fe_20-Co_20-Ni_20-Cr_20-Mn_20/sample_0001/openmx_slab/input.scfout
```

## 8. Extract Surface d-Orbital Hamiltonian

```powershell
python extract_openmx_hamiltonian.py --scfout dataset/surfaces/Fe_20-Co_20-Ni_20-Cr_20-Mn_20/sample_0001/openmx_slab/input.scfout --dat dataset/surfaces/Fe_20-Co_20-Ni_20-Cr_20-Mn_20/sample_0001/openmx_slab/input.dat --top-atoms dataset/surfaces/Fe_20-Co_20-Ni_20-Cr_20-Mn_20/sample_0001/metadata/top_atoms.jsonl --output dataset/surfaces/Fe_20-Co_20-Ni_20-Cr_20-Mn_20/sample_0001/openmx_slab/hamiltonian_d_surface.npz --dataset-root dataset --surface-id Fe_20-Co_20-Ni_20-Cr_20-Mn_20 --sample-id sample_0001
```

## Final Sample Contents

```text
sample_0001/
  manifest.json
  structures/
    00_initial_sqs.cif
    01_relaxed_slab.cif
  metadata/
    top_atoms.jsonl
    atom_grid.npy
    fcc_sites.jsonl
    site_grid.npy
  adsorbates/
    N/
      site_0001/
        00_initial_adsorbate.cif
        01_relaxed_adsorbate.cif
        adsorption_energy.json
      ...
  openmx_slab/
    input.dat
    input.scfout
    hamiltonian_d_surface.npz
    hamiltonian_d_surface.npz.basis.jsonl
```
