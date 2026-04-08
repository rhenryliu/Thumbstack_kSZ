# Catalogue Pipeline

This directory contains scripts and configs for building DESI Y3 galaxy catalogues from raw
FITS files. The outputs are consumed by two separate pipelines:

- **kSZ stacking** — this repo (`ThumbStack_kSZ`), via ThumbStack-format `.txt` files
- **Weak lensing** — `~/projects/DESIxHSC-Lensing`, via lensing-format `.fits` files

## Overview

Building a new catalogue requires three steps across two Python environments:

| Step | Script(s) | Environment |
|------|-----------|-------------|
| 1. Raw DESI FITS → intermediate CSV | `prepare_catalogue_yaml.py` or `prepare_catalogue_all_yaml.py` | `cosmodesi_dr1` venv |
| 2. CSV → ThumbStack `.txt` | `convert_to_TS_catalog_allbins.py` + `convert_to_TS_catalog.py` | `rhliu_tSZ` conda env |
| 3. CSV → lensing `.fits` | `~/projects/DESIxHSC-Lensing/data_processing/convert_sz_to_lensing_cat.py` | `cosmodesi_dr1` venv |

Steps 1 and 3 share the same environment. Step 2 must be separate because ThumbStack's
dependencies are incompatible with `cosmodesi`.

---

## Step 1 — Process raw DESI catalogues

Two scripts handle different input formats:

- **`prepare_catalogue_yaml.py`** — For catalogues taken directly from the DESI survey,
  distributed as separate NGC and SGC FITS files. Computes LOS peculiar velocities from the
  displacement between pre- and post-BAO-reconstruction positions using `cosmoprimo`. This is
  the standard script for LRG and BGS catalogues from the official DESI data releases.
- **`prepare_catalogue_all_yaml.py`** — For catalogues that have already been processed
  externally and are provided as a single FITS file with a reconstructed LOS velocity column
  `vR` already present. No velocity computation is performed and there is no NGC/SGC split.
  Currently used for Boryana Hadzhiyska's BGS CIGALE mass catalogue, which combines both
  caps and includes pre-computed velocities.

Run from `desi/spec_Y3/catalogue/`:

```bash
source /global/common/software/desi/users/adematti/cosmodesi_environment.sh dr1

python -u prepare_catalogue_yaml.py -p ./configs/prepare_cat_LRG.yaml
```

### Output files

Saved to `save_dir` (from config):

- `{field}_catalog_Y3_{filter_type}.csv` — full-z catalogue (`field` = `full`, `NGC`, `SGC`)
- `{field}_zbin{N}_{filter_type}.csv` — per-redshift-bin catalogue

`prepare_catalogue_all_yaml.py` only produces `full` field outputs (no NGC/SGC split).

### Config parameters

#### `source` section — for `prepare_catalogue_yaml.py`

| Parameter | Type | Description |
|-----------|------|-------------|
| `main_directory` | str | Root path of the DESI DA2 catalogue tree |
| `pre_rec_dir` | str | Relative path (from `main_directory`) to pre-reconstruction FITS files |
| `post_rec_dir` | str | Relative path to post-reconstruction FITS files |
| `NGC_fn` | str | Filename of the NGC clustering catalogue |
| `SGC_fn` | str | Filename of the SGC clustering catalogue |

#### `source` section — for `prepare_catalogue_all_yaml.py`

| Parameter | Type | Description |
|-----------|------|-------------|
| `cat_fn` | str | Full path to the single input FITS file (must contain a `vR` column) |

#### `processing` section

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `cat_type` | str | `'LRG'` | Galaxy type: `'LRG'` or `'BGS'`. Used in output file naming. |
| `filter_type` | str | `'nopairs'` | Sky/quality filter applied to the catalogue (see below) |
| `physical_sep` | float | `1.0` | Minimum physical separation in **Mpc** for close-pair removal. Only used when `filter_type: nopairs`. |
| `zeroV` | bool | `false` | If true, subtracts the mean LOS peculiar velocity separately for each field (NGC/SGC). Only supported in `prepare_catalogue_yaml.py`. |
| `zbins` | list of `[z_min, z_max]` | — | Redshift bin edges. Each pair produces a separate output CSV. |

**Filter type options:**

| Value | ACT mask | Pair pruning | Notes |
|-------|----------|--------------|-------|
| `unfiltered` | No | No | All galaxies with valid positions |
| `no_src_with_cluster_mask` | Yes | No | Removes objects outside ACT footprint or near point sources/clusters |
| `nopairs` | Yes | Yes | ACT mask + removes close pairs within `physical_sep` Mpc (greedy, catalog-order priority) |

#### `save` section

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `save_dir` | str | — | Output directory for CSV files |
| `save_format` | str | `'csv'` | Output format: `'csv'` or `'txt'` |

### Existing configs

| Config | Script | `cat_type` | `filter_type` | Notes |
|--------|--------|------------|---------------|-------|
| `prepare_cat_LRG.yaml` | `prepare_catalogue_yaml.py` | LRG | `no_src_with_cluster_mask` | Primary LRG catalogue |
| `prepare_cat_LRG_zeroV.yaml` | `prepare_catalogue_yaml.py` | LRG | `no_src_with_cluster_mask` | LRG with mean velocity zeroed (null variant) |
| `prepare_cat_LRG_null.yaml` | `prepare_catalogue_yaml.py` | LRG | `nopairs` | Random catalogue for null tests (`NGC_fn`/`SGC_fn` point to `.ran.fits` files) |
| `prepare_cat_BGS.yaml` | `prepare_catalogue_yaml.py` | BGS | `nopairs` | Primary BGS catalogue |
| `prepare_cat_BGS_cigale.yaml` | `prepare_catalogue_all_yaml.py` | BGS | `nopairs` | BGS from Boryana's CIGALE mass catalogue (single-file input) |
| `prepare_cat_BGS_null.yaml` | `prepare_catalogue_yaml.py` | BGS | `nopairs` | Random catalogue for null tests |

`prepare_cat_LRG_2.yaml`, `prepare_cat_LRG_3.yaml`, `prepare_cat_BGS_2.yaml`, and
`prepare_cat_BGS_3.yaml` are older configs retained for reference; their filter variants are
now configurable via CLI or a fresh config. `prepare_catalogue.py` and
`prepare_catalogue_nopairs.py` are deprecated hardcoded predecessors — do not use.

---

## Step 2 — Convert to ThumbStack format

Use the same config file as step 1, but activate `rhliu_tSZ`:

```bash
module load python
conda activate rhliu_tSZ

python -u convert_to_TS_catalog_allbins.py -p ./configs/prepare_cat_LRG.yaml
python -u convert_to_TS_catalog.py -p ./configs/prepare_cat_LRG.yaml
```

`convert_to_TS_catalog_allbins.py` converts the full-z CSVs (`full_catalog_Y3_*.csv`,
`NGC_catalog_Y3_*.csv`, `SGC_catalog_Y3_*.csv`) into ThumbStack `.txt` files. \
`convert_to_TS_catalog.py` converts the per-z-bin CSVs (`{field}_zbin{N}_*.csv`).

Output files (written to the same `save_dir`):

- `DESIY3_{cat_type}/catalog_{field}_{filter_type}.txt` — full-z catalogue
- `DESIY3_{cat_type}_z{N}/catalog_{field}_{filter_type}.txt` — per-z-bin

These `.txt` files are what `make_stacks_yaml.py` reads.

---

## Step 3 — Convert to lensing format

Edit and run `~/projects/DESIxHSC-Lensing/data_processing/convert_sz_to_lensing_cat.py`
(cosmodesi_dr1 environment). See the `data_processing/README.md` in that repo for details.

---

## SLURM scripts

| Script | Purpose |
|--------|---------|
| `make_catalogues.sh` | Steps 1+2 for LRG. Step 1 lines are commented out; uncomment as needed. |
| `make_catalogues2.sh` | Step 1 for BGS or null configs. Uses cosmodesi_dr1. |
| `make_catalogues_cigale.sh` | Step 2 for CIGALE BGS configs. Uses rhliu_tSZ. |

All scripts use `--account=m3058` and `--qos=debug`. Upgrade to `regular` for large runs.

---

## Output location

Catalogues land in subdirectories of `/pscratch/sd/r/rhliu/projects/Weak_lensing/desi/spec_Y3/`
as configured in each config's `save_dir`. Typical layout:

```
spec_Y3/
├── LRG_catalogues/
│   ├── full_catalog_Y3_no_src_with_cluster_mask.csv
│   ├── full_zbin1_no_src_with_cluster_mask.csv
│   ├── ...
│   ├── DESIY3_LRG/catalog_full_no_src_with_cluster_mask.txt
│   └── DESIY3_LRG_z1/catalog_full_no_src_with_cluster_mask.txt
├── BGS_catalogues/
└── BGS_CIGALE/
```
