# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Overview

This repo contains a DESI-specific analysis pipeline for measuring the kinematic Sunyaev-Zel'dovich (kSZ) effect by stacking CMB temperature maps at galaxy positions. It wraps the core **ThumbStack** library (a separate repo at `/global/homes/r/rhliu/projects/repos/ThumbStack`) and adds DESI catalogue preparation, filtering, and stacking scripts.

## kSZ filter terminology for the DESIxHSC-Lensing f_gas paper (2026-10-23)

- **v2 filter** = option C = `DSigma_lensing`, the lensing-matched kSZ filter (the lensing-bin
  average of ΔΣ on annuli matched to the lensing bins); in SimulationStacker it includes filtering the
  simulations on ThumbStack-like 0.25′ stamps, as here for the data. The production filter from
  2026-10-23.
- **v1 filter** = `DSigma`, the published filter (mean Σ inside R minus a 0.75′ ring outside R).
- Write "the v2 filter (option C, `DSigma_lensing`)" or "the v1 filter (`DSigma`)", with the code
  name, in docs, reports and captions. Always write "lensing v2.1" and "the v2 filter" in full, never
  a bare "v2", since the lensing also has a v2. Never write "paper filter", "the paper's filter",
  or the interim names "default filter" and "old filter".
- Option A (the v1 filter divided by the simulations' B′ beam factor) was rejected on 2026-10-23:
  the simulation spread in B′ is too large for a model-independent correction. Do not use or
  suggest it.
- Make or show v1-filter results only when asked for that comparison, and keep them in folders
  whose names say the filter, never next to v2-filter outputs.
- Here the filter is `stack.filterType` (`'DSigma_lensing'` or `'DSigma'`; `diskring` is the CAP
  filter); the v2 filter's configs have `dsigma_lensing` in their names.

## f_gas version terminology for the DESIxHSC-Lensing f_gas paper (2026-10-25)

The versions name which probes carry the DESI completeness weights; the f_gas ratios are formed in
DESIxHSC-Lensing. kSZ paths are under `/pscratch/sd/r/rhliu/projects/Weak_lensing/ksz_measurements/ACTxDESI/`.

- **v2.1** = neither the kSZ nor the lensing is weighted: the unweighted v2-filter kSZ stacks in
  `spec_Y3_DSigma_lensing/` (estimator `ksz_uniformweight`), over lensing v2.1. The default.
- **v2.2** = both are weighted: the completeness-weighted kSZ (estimator `ksz_compweight`, configs in
  `desi/spec_Y3/configs/compweight/`) in `spec_Y3_DSigma_lensing_compweight/`, over lensing v2.2.
  That folder also holds an unweighted restack of the same objects; v2.1 reads
  `spec_Y3_DSigma_lensing/`. Whether v2.2 replaces v2.1 is not decided.
- Weighting only one probe is wrong: both or neither. The lensing-only weighting (lensing v2.2 over
  the unweighted kSZ) is kept in DESIxHSC-Lensing only for a possible appendix, tagged
  `mixed_weights`. Never call it v2.2; write "lensing-only completeness weighting (`mixed_weights`)".
- Details: `docs/completeness_weighted_ksz.md` (plan) and `docs/completeness_weighted_ksz_results.md`
  (results and null tests) here; DESIxHSC-Lensing `docs/completeness_weights.md` and
  `docs/lensing_v2.2.md`.

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

Other overrides: `--cat-type` (LRG/BGS) and `--cat-names` (space-separated list).

SLURM wrappers: `stack_full.sh`, `stack_NGC.sh`, `stack_SGC.sh`, `stack_CAP.sh`, `stack_dsigma.sh`, etc.

## Testing

pytest suite, run in the `rhliu_tSZ` environment from the repository root:

```bash
python -m pytest tests/
```

- `tests/test_dsigma_lensing.py`: Gaussian known answers for the `DSigma_lensing` filter
  (`dsigma_lensing.py` and its branches in `thumbstack.py`).
- `tests/test_ksz_compweight.py`: the `ksz_compweight` estimator against an independent reference,
  and `desi/spec_Y3/weight_utils.py`.

`desi/spec_Y3/checks/` holds one-off checks kept as records (most behind
`docs/completeness_weighted_ksz_results.md`; `data_jackknife.py` is the v2-filter data-stack
jackknife); each docstring names its recorded run. `desi/spec_Y3/compweight_nulls.py`
evaluates the completeness-weighted null tests from the outputs of `configs/compweight/`.

## Architecture

### ThumbStack dependency

Core classes (`ThumbStack`, `Catalog`, `cmbMap`, `UnivMariana`, `MassConversionKravtsov14`, etc.) are imported from the upstream repo at `/global/homes/r/rhliu/projects/repos/ThumbStack` (or `/global/homes/r/rhliu/projects/repos/ThumbStack_kSZ` for local copies). All scripts add one of these paths via `sys.path.append`. The classes are imported with `import *` (intentional, preserved from upstream).

### Config-driven stacking

`make_stacks_yaml.py` reads a YAML with these top-level sections:
- `stack`: `filterType`, `save`, `zeroV`; optional aperture grid `rApMinArcmin`, `rApMaxArcmin`, `nRAp`
  (without them ThumbStack keeps its hardcoded radii); optional `shuffle`, `shuffle_seed`
  - `filterType` values in `thumbstack.py`: `DSigma_lensing` (the v2 filter), `DSigma` (the v1
    filter), `diskring` (CAP), `disk`, `ring`, `meanring`, `upsilon`
  - completeness-weighted stacks (all optional; without them the run is unchanged): `estimators`
    (e.g. `['ksz_uniformweight', 'ksz_compweight']`), `filtering_dir` (read-only per-object filtering
    of a production stack), `weights` (columns, see `weight_utils.load_weights`), `vshuffle` and
    `n_vshuffle`, `save_samples`. `vshuffle` needs `estimators`, and `weights` cannot be combined
    with the `mass` section; both raise errors
- `data`: paths, `cat_type` (LRG/BGS), `field` (NGC/SGC/full), `filter_cut` (unfiltered/no_src_with_cluster_mask/nopairs), `cat_names` list
- `extra`: whether to save individual per-object stacks and where
- `mass` (optional): stellar-mass cuts or bins applied at stacking time (`mass_filter_utils.py`),
  e.g. `col`, `strategy`, `n_bins`, `bin_idx`, `drop_zero_mass`

YAML configs for stacking are in `desi/spec_Y3/configs/`. Catalogue prep configs are in `desi/spec_Y3/catalogue/configs/`.
Subfolders of `desi/spec_Y3/configs/`:
- `compweight/`: completeness-weighted stacks and their nulls, for the v2 filter (`*_dsigma_lensing_*`)
  and the v1 filter (`*_dsigma_*`, outputs in `spec_Y3_DSigma_compweight*/` folders); `checks/`:
  known-answer and no-new-keys (`defaults_*`) configs
- `masscut/`: stellar-mass bins and cuts
- `null/`: random-position and shuffled-velocity nulls
- `ext10/`, `DES/`: extended-aperture and DES/non-DES footprint variants
- `tests/`: the median-estimator diagnostic (`make_stacks_median.py`); `old/`: superseded configs

### Catalogue preparation config

`prepare_cat_LRG.yaml` (and BGS equivalents) controls:
- `source`: DESI DR file paths for pre-/post-reconstruction catalogues
- `processing`: `cat_type`, `filter_type`, redshift bins `zbins`
- `save`: output directory and format

Filter types for catalogue prep:
- `nopairs`: remove close pairs within a physical separation (default 0.5 Mpc)
- `no_src_with_cluster_mask`: source/cluster masking
- `unfiltered`: no additional cuts

### Key modules in this repo

- `thumbstack.py`: this repo's fork of the `ThumbStack` class; the stacking filters (including
  `DSigma_lensing`) and the `ksz_compweight` estimator are implemented here
- `dsigma_lensing.py`: weights of the `DSigma_lensing` (v2) filter
- `universe.py`: cosmological utilities (`UnivMariana`)
- `mass_conversion.py`: M*–Mh relations
- `cmb.py`, `cmbMap.py`: CMB map loading and handling
- `flat_map.py`: flat-sky map operations
- `desi/spec_Y3/catalogue/remove_close_pairs.py`: pair-pruning logic
- `desi/spec_Y3/catalogue/catalog_utils.py`: shared catalogue utilities
- `desi/spec_Y3/weight_utils.py`: per-object completeness weights for `ksz_compweight`
- `desi/spec_Y3/mass_filter_utils.py`: stellar-mass filters applied at stacking time (a partial copy
  lives in DESIxHSC-Lensing `src/mass_filter_utils.py`)

### Output structure

Stacking outputs land in subdirectories of the `output_dir` specified in the config, organised by catalogue name and CMB map name. The `save=True` flag in the stacking config recomputes and saves temperature decrements; `save=False` loads previously saved ones.

## Pre-push code review hook

`.claude/hooks/pre-push` runs an autonomous Claude code review before pushes. It is opt-in:
```bash
REVIEW=1 git push
```
Without `REVIEW=1`, the hook exits immediately without reviewing.
