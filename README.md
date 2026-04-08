# ThumbStack — DESI Y3 kSZ Stacking Pipeline

ThumbStack produces stacked maps and radial profiles from catalogs of object positions and
CMB maps. It is designed for thermal and kinematic Sunyaev-Zel'dovich (tSZ/kSZ)
measurements, and outputs 2D stacked maps and profiles for different spatial filters (e.g.,
aperture photometry, disk-minus-ring), along with covariance estimates via bootstrap.

This repo is a fork of the upstream ThumbStack library with DESI-specific catalogue
preparation and stacking scripts added in `desi/spec_Y3/`.

ThumbStack was originally used in the following publication:
https://ui.adsabs.harvard.edu/abs/2020arXiv200905557S/abstract

Please feel free to use this code in your work, and contact me with any questions or
suggestions (eschaan@lbl.gov). If you use it in your publication, please cite the paper above.

---

## DESI Y3 kSZ Pipeline

Measures the kSZ effect by stacking ACT DR6 CMB temperature maps at DESI Y3 galaxy
positions, weighted by reconstructed LOS peculiar velocities from BAO reconstruction.

### Data

- **CMB maps**: ACT DR6 harmonic-ILC temperature maps and masks from
  `/pscratch/sd/b/boryanah/ACTxDESI/ACT/`
- **Galaxy catalogues**: DESI Y3 LRGs and BGS (spectroscopic), processed into
  ThumbStack-format `.txt` files under
  `/pscratch/sd/r/rhliu/projects/Weak_lensing/desi/spec_Y3/`

### Environment

```bash
module load python
conda activate rhliu_tSZ
```

Catalogue preparation (step 1 of the pipeline) requires a separate `cosmodesi_dr1`
environment — see [`desi/spec_Y3/catalogue/README.md`](desi/spec_Y3/catalogue/README.md).

### Pipeline overview

**Stage 1 — Catalogue preparation** (`desi/spec_Y3/catalogue/`)

Converts raw DESI FITS files into ThumbStack-format `.txt` catalogues, computing LOS
peculiar velocities from BAO reconstruction displacements along the way. See
[`desi/spec_Y3/catalogue/README.md`](desi/spec_Y3/catalogue/README.md) for the full
step-by-step workflow, config parameter reference, and SLURM scripts.

**Stage 2 — Stacking** (`desi/spec_Y3/`)

Run with:

```bash
python -u make_stacks_yaml.py -p ./configs/LRG_diskring_full.yaml

# Override config values at the command line (no need to copy the YAML for quick tests):
python -u make_stacks_yaml.py -p ./configs/LRG_diskring_full.yaml \
    --filterType diskring --field NGC --filter-cut nopairs
```

`make_stacks_yaml.py` loads the galaxy catalogue and ACT map, then runs `ThumbStack` to
compute the mean CMB temperature decrement as a function of angular scale (the kSZ profile).
Configs live in `desi/spec_Y3/configs/`. Key config sections:

- **`stack`**: `filterType` (`diskring` or `DSigma`), `save` (recompute vs load cached
  decrements), `zeroV` (subtract mean velocity at runtime), `shuffle` (null test via velocity
  shuffling)
- **`data`**: `cat_type` (LRG/BGS), `field` (NGC/SGC/full), `filter_cut`, `cat_names`
  (list of ThumbStack catalogue names to load), input/output paths
- **`extra`**: optional per-object stack saving

SLURM wrappers: `stack_full.sh`, `stack_NGC.sh`, `stack_SGC.sh`, `stack_dsigma.sh`,
`stack_CAP.sh` — all use `--account=desi`.

### Output location

```
/pscratch/sd/r/rhliu/projects/Weak_lensing/ksz_measurements/ACTxDESI/spec_Y3/
└── {extra_str}/
    └── {cat_name}_{cmbMap_name}/
```

### Related repo

`~/projects/DESIxHSC-Lensing` measures weak lensing (ΔΣ) with HSC using the same DESI
catalogues, allowing a joint kSZ + lensing analysis.

### Further reading

See [`CLAUDE.md`](CLAUDE.md) for architecture details and module descriptions.
