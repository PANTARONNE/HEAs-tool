# Minimal Input Workflow

This is the shortest current workflow for the example composition
`Co_13-Cr_13-Fe_13-Mn_12-Ni_13`.

## Required Manual Inputs

```text
elements              Fe Co Ni Cr Mn
ratios (optional)     1 1 1 1 1        # random when omitted
relaxed_slab_cif      test_sqs-opt.cif
adsorbate             N
openmx_data_path      DFT_DATA19
openmx_scfout         input.scfout
adsorption_energy_eV  one value per site
```

The composition ID is derived automatically from the generated atom counts. It
is not a command-line input to `create-sample`.

## 1. Generate and Register the Structure

```powershell
python hea_dataset.py create-sample --root dataset --elements Fe Co Ni Cr Mn --ratios 1 1 1 1 1
```

Created paths (only the initial SQS structure at this stage):

```text
dataset/index.sqlite
dataset/dataset_manifest.json
dataset/Co_13-Cr_13-Fe_13-Mn_12-Ni_13/manifest.json
dataset/Co_13-Cr_13-Fe_13-Mn_12-Ni_13/structures/00_initial_sqs.cif
```

## 1b. Register the Relaxed Structure

After external relaxation of `00_initial_sqs.cif`:

```powershell
python hea_dataset.py record-relaxed --root dataset --surface-id Co_13-Cr_13-Fe_13-Mn_12-Ni_13 --relaxed-cif test_sqs-opt.cif
```

Created path:

```text
dataset/Co_13-Cr_13-Fe_13-Mn_12-Ni_13/structures/01_relaxed_slab.cif
```

## 2. Index Surface Atoms

```powershell
python hea_dataset.py index-surface --root dataset --surface-id Co_13-Cr_13-Fe_13-Mn_12-Ni_13
```

Created paths:

```text
dataset/Co_13-Cr_13-Fe_13-Mn_12-Ni_13/metadata/top_atoms.jsonl
dataset/Co_13-Cr_13-Fe_13-Mn_12-Ni_13/metadata/atom_grid.npy
```

## 3. Detect FCC Sites

```powershell
python hea_dataset.py detect-sites --root dataset --surface-id Co_13-Cr_13-Fe_13-Mn_12-Ni_13
```

Created paths:

```text
dataset/Co_13-Cr_13-Fe_13-Mn_12-Ni_13/metadata/fcc_sites.jsonl
dataset/Co_13-Cr_13-Fe_13-Mn_12-Ni_13/metadata/site_grid.npy
```

## 4. Create Initial Adsorbate Structures

```powershell
python hea_dataset.py create-adsorbate-records --root dataset --surface-id Co_13-Cr_13-Fe_13-Mn_12-Ni_13 --adsorbate N
```

The command uses the registered `01_relaxed_slab.cif` and creates one directory
per detected site:

```text
dataset/Co_13-Cr_13-Fe_13-Mn_12-Ni_13/adsorbates/N/site_0001/00_initial_adsorbate.cif
dataset/Co_13-Cr_13-Fe_13-Mn_12-Ni_13/adsorbates/N/site_0001/adsorption_energy.json
dataset/Co_13-Cr_13-Fe_13-Mn_12-Ni_13/adsorbates/N/site_0002/00_initial_adsorbate.cif
dataset/Co_13-Cr_13-Fe_13-Mn_12-Ni_13/adsorbates/N/site_0002/adsorption_energy.json
...
```

After external relaxation, save each relaxed structure as:

```text
dataset/Co_13-Cr_13-Fe_13-Mn_12-Ni_13/adsorbates/N/site_0001/01_relaxed_adsorbate.cif
```

## 5. Record One Adsorption Energy

```powershell
python hea_dataset.py record-energy --root dataset --surface-id Co_13-Cr_13-Fe_13-Mn_12-Ni_13 --adsorbate N --site-id 1 --relaxed-cif N-site-0001-relaxed.cif --energy -1.23
```

The command validates the relaxed structure and copies it to
`adsorbates/N/site_0001/01_relaxed_adsorbate.cif`. Adsorbate displacement must
not exceed 2.0 Å by default. Use `--max-adsorbate-displacement` to change the
threshold. Repeat for every site.

## 6. Generate the Clean-Slab OpenMX Input

```powershell
python cif_to_openmx.py dataset/Co_13-Cr_13-Fe_13-Mn_12-Ni_13/structures/01_relaxed_slab.cif --data-path DFT_DATA19 -o dataset/Co_13-Cr_13-Fe_13-Mn_12-Ni_13/openmx_slab/input.dat
```

Place the OpenMX output at:

```text
dataset/Co_13-Cr_13-Fe_13-Mn_12-Ni_13/openmx_slab/input.scfout
```

## 7. Extract the Surface Hamiltonian

```powershell
python hea_dataset.py extract-hamiltonian --root dataset --surface-id Co_13-Cr_13-Fe_13-Mn_12-Ni_13 --dat dataset/Co_13-Cr_13-Fe_13-Mn_12-Ni_13/openmx_slab/input.dat --scfout dataset/Co_13-Cr_13-Fe_13-Mn_12-Ni_13/openmx_slab/input.scfout
```

`top_atoms.jsonl` and the default output path are resolved automatically from
the surface ID. `--dat` and `--scfout` must be supplied explicitly.

## Final Directory Structure

```text
dataset/
  index.sqlite
  dataset_manifest.json
  Co_13-Cr_13-Fe_13-Mn_12-Ni_13/
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
