# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Overview

This repo contains a DESI-specific analysis pipeline for measuring the kinematic Sunyaev-Zel'dovich (kSZ) effect by stacking CMB temperature maps at galaxy positions. It wraps the core **ThumbStack** library (a separate repo at `/global/homes/r/rhliu/projects/repos/ThumbStack`) and adds DESI catalogue preparation, filtering, and stacking scripts.

## Environment

- Conda env: `rhliu_tSZ` (activate before running any scripts)
- SLURM accounts: `desi` (default), `m3058` (catalogue prep jobs)
- CMB map data: `/pscratch/sd/b/boryanah/ACTxDESI/ACT/`
- Output data: `/pscratch/sd/r/rhliu/projects/Weak_lensing/ksz_measurements/ACTxDESI/spec_Y3/`
- Catalogue data: `/pscratch/sd/r/rhliu/projects/Weak_lensing/desi/spec_Y3/`

## Running the pipeline

The active pipeline lives under `desi/spec_Y3/`. There are two stages:

### 1. Catalogue preparation (run from `desi/spec_Y3/catalogue/`)

```bash
# Prepare intermediate catalogues (requires cosmodesi environment)
python -u prepare_catalogue_yaml.py -p ./configs/prepare_cat_LRG.yaml

# Convert to ThumbStack-format .txt and per-z-bin .csv files
python -u convert_to_TS_catalog.py -p ./configs/prepare_cat_LRG.yaml
python -u convert_to_TS_catalog_allbins.py -p ./configs/prepare_cat_LRG.yaml
```

SLURM wrapper: `desi/spec_Y3/catalogue/make_catalogues.sh`

### 2. Stacking (run from `desi/spec_Y3/`)

The preferred entry point is `make_stacks_yaml.py` (supports CLI overrides); `dsigma_thumbstack_yaml.py` is the older equivalent.

```bash
python -u make_stacks_yaml.py -p ./configs/LRG_diskring_full.yaml
python -u make_stacks_yaml.py -p ./configs/LRG_dsigma_full.yaml

# CLI overrides (skip editing the YAML for quick tests):
python -u make_stacks_yaml.py -p ./configs/LRG_diskring_full.yaml \
    --filterType diskring --field NGC --filter-cut unfiltered
```

SLURM wrappers: `stack_full.sh`, `stack_NGC.sh`, `stack_SGC.sh`, `stack_CAP.sh`, `stack_dsigma.sh`, etc.

## Architecture

### ThumbStack dependency

Core classes (`ThumbStack`, `Catalog`, `cmbMap`, `UnivMariana`, `MassConversionKravtsov14`, etc.) are imported from the upstream repo at `/global/homes/r/rhliu/projects/repos/ThumbStack` (or `/global/homes/r/rhliu/projects/repos/ThumbStack_kSZ` for local copies). All scripts add one of these paths via `sys.path.append`. The classes are imported with `import *` (intentional, preserved from upstream).

### Config-driven stacking

`make_stacks_yaml.py` reads a YAML with three top-level sections:
- `stack`: `filterType` (`diskring` or `DSigma`), `save`, `zeroV`
- `data`: paths, `cat_type` (LRG/BGS), `field` (NGC/SGC/full), `filter_cut` (unfiltered/no_src_with_cluster_mask/nopairs), `cat_names` list
- `extra`: whether to save individual per-object stacks and where

YAML configs for stacking are in `desi/spec_Y3/configs/`. Catalogue prep configs are in `desi/spec_Y3/catalogue/configs/`.

### Catalogue preparation config

`prepare_cat_LRG.yaml` (and BGS equivalents) controls:
- `source`: DESI DR file paths for pre-/post-reconstruction catalogues
- `processing`: `cat_type`, `filter_type`, redshift bins `zbins`
- `save`: output directory and format

Filter types for catalogue prep:
- `nopairs`: remove close pairs within a physical separation (default 0.5 Mpc/h)
- `no_src_with_cluster_mask`: source/cluster masking
- `unfiltered`: no additional cuts

### Key modules in this repo

- `universe.py`: cosmological utilities (`UnivMariana`)
- `mass_conversion.py`: M*–Mh relations
- `cmb.py`, `cmbMap.py`: CMB map loading and handling
- `flat_map.py`: flat-sky map operations
- `desi/spec_Y3/catalogue/remove_close_pairs.py`: pair-pruning logic
- `desi/spec_Y3/catalogue/catalog_utils.py`: shared catalogue utilities

### Output structure

Stacking outputs land in subdirectories of the `output_dir` specified in the config, organised by catalogue name and CMB map name. The `save=True` flag in the stacking config recomputes and saves temperature decrements; `save=False` loads previously saved ones.
