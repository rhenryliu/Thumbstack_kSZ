# Completeness-weighted kSZ stacks: context and plan

*Written 2026-10-23 by the session that made lensing v2.2 in DESIxHSC-Lensing, as the brief for
the session that will do this work. Status: **plan approved, not started.** Nothing in this repo
has been changed for it yet. Whether the paper adopts completeness weighting is **not decided**:
this work produces the weighted kSZ stacks and their null tests so that the choice can be made.
The full literature review and the measurements behind this plan are in
`~/projects/DESIxHSC-Lensing-v2.2/docs/completeness_weights.md` (branch `lensing-v2.2`; on `main`
once merged) — read it first.*

---

## 1. The project and why this is needed

**The measurement.** The gas fraction f_gas(θ) = ΔΣ_kSZ(θ) / ΔΣ_lens(θ) (Liu et al. 2026,
arXiv:2609.22466) divides a velocity-weighted kSZ stack of ACT DR6 maps on DESI DR2 galaxies, made
in this repo, by the galaxy–galaxy lensing ΔΣ of the same galaxy samples measured with HSC Y3
sources in DESIxHSC-Lensing (`dsigma`). The ratio is formed in DESIxHSC-Lensing
(`scripts/plot_scripts/compare_combined_{self,BGS}.py`), which reads this repo's
`{filterType}_ksz_uniformweight_measured.txt` and `cov_{filterType}_ksz_uniformweight_bootstrap.txt`
and divides by the velocity reconstruction coefficient r_V (`r_val`: 0.65 for LRG from Qu et al.
2026, arXiv:2604.19744; 0.64 for BGS from Hadzhiyska et al. 2026, arXiv:2604.19745). The data are
compared with hydrodynamical simulations in SimulationStacker.

**Samples.** LRG bins 1–4 (0.4–0.6, 0.6–0.8, 0.8–0.95, 0.95–1.1) and BGS (one bin, CIGALE stellar
masses). The kSZ samples are ACT-masked (`no_src_with_cluster_mask`). For BGS the kSZ and lensing
use the same 2,224,546-row CIGALE file. For LRG the lensing uses the full 1.05M-lens bin-1 sample
and the kSZ its 597k ACT-masked subset, by design.

**Filters.** `DSigma` (the paper's) and `DSigma_lensing` (lensing-matched annuli, "option C"; the
current integration). Both have per-object products saved for every stack used below.

**Lensing versions** (DESIxHSC-Lensing docs `lensing_v2*.md`):
- v2: HSC shear calibration (m, m_sel), n(z) from Choppin de Janvry et al., bug fixes.
- **v2.1: v2 with the lens-magnification term subtracted; the current baseline.**
- **v2.2: v2.1 with DESI completeness weights on the lenses** (`WEIGHT_COMP × WEIGHT_ZFAIL` on
  top of `WEIGHT_SYS`) and `FRAC_TLOBS_TILES` on the randoms. This is the DESI DR1 lensing
  convention (Heydenreich et al. 2025, arXiv:2506.21677). It raises ΔΣ by ×1.035–1.10 (LRG bin 1)
  and ×1.14 (1′) to ×1.49 (6′) for BGS. A test; it does not replace v2.1.

**Why the kSZ must be weighted too.** DESI misses targets where they crowd (BGS completeness is
below 10 % in cluster cores, Smith et al. 2019), so the observed sample under-represents dense
environments and satellites; the completeness weight restores the full target sample. In both
probes a satellite carries its host's halo, so kSZ and ΔΣ must have the same satellite mix:
**weight both or neither.** v2.2 weights only the lensing, so its f_gas is not to be used.
B. Hadzhiyska (collaborator) argues for weighting both; R. H. Liu's concern is that uniform kSZ
weighting is the field's standard, and that extra weights break the zero mean velocity and could
let tSZ or dust leak in.

## 2. What is already known

**The weights.** DA2 `loa-v1/LSScats/v1.1/nonKP` (where both lens samples come from):
`WEIGHT_COMP` = 1/`FRACZ_TILELOCID` exactly, integer valued, 1 for 82 % (BGS) / 88 % (LRG) of
galaxies. `LSScats/v2/PIP`: `WEIGHT_COMP` = 129/(1 + k), k the number of 128 alternative
fibre-assignment runs with a fibre (`BITWEIGHTS`): the individual inverse probability. The two
correlate per galaxy at 0.48 (LRG) and 0.62 (BGS).

**Literature, in brief.** DESI lensing weights its lenses as standard (Lange et al. 2024,
arXiv:2404.09397, recommend it for all galaxy–galaxy lensing: unweighted BGS ΔΣ is up to 30 % low,
LRG ~5 %). No velocity-weighted kSZ stack in the literature applies completeness weights,
including the DESI DR2 kSZ papers, which name fibre undersampling as an uncorrected caveat. The
only DESI kSZ analysis with them is Chaussidon et al. 2026 (arXiv:2604.04867), a power spectrum.
Satellites are 20–30 % of DESI-like samples but carry 50–80 % of the stacked kSZ signal
(Ondaro-Mallea et al. 2026, arXiv:2610.17200).

**A weighted kSZ check already exists** (DESIxHSC-Lensing-v2.2,
`scratch/v2.2_ksz_weights/ksz_weight_check.py` → `results_ksz_weights.txt`, job 60160843). It
rebuilds the estimator from this repo's saved per-object arrays with `src/ksz_reconstruct.py` (the
unweighted stack reproduces the saved profiles to 3 × 10⁻¹⁴). Weighted estimator:

  T̂_w(R) = σ_w / r_V · Σ w_i v′_i T_i(R) / Σ w_i v′_i²,  v′_i = v_i − Σ w v / Σ w,
  σ_w² = Σ w (v − v̄_w)² / Σ w,

over the objects that survive the catalogue mask, with w = `WEIGHT_COMP × WEIGHT_ZFAIL`.
Results (v1.1 weights, weighted/unweighted, paired bootstrap):

| Sample, filter | kSZ | f_gas both weighted / both unweighted |
|---|---|---|
| BGS, `DSigma_lensing` | 1.22 ± 0.08 | 0.95 ± 0.07 |
| BGS, `DSigma` | 1.23 ± 0.08 | 1.03 ± 0.07 |
| LRG bin 1, `DSigma_lensing` | 0.99 ± 0.04 | 0.94 ± 0.04 |
| LRG bin 1, `DSigma` | 1.00 ± 0.04 | 0.96 ± 0.04 |

- PIP weights give the same kSZ response (BGS 1.21 / 1.22, LRG bin 1 0.98 / 0.99).
- **The leakage concern:** the weighted mean velocity is −0.0035 σ_v (BGS) and +0.0127 σ_v
  (LRG bin 1). Not re-centring on it leaks ≤ 0.24σ per aperture; re-centring removes it exactly.
  200 weighted shuffles: mean ≤ 0.15σ, as unweighted. Weighted null scatter ×1.13–1.26 (BGS),
  ×1.02–1.11 (LRG bin 1).
- **Side finding:** `make_stacks_yaml.py` zeroes `vR` before the catalogue mask, so the stacked
  objects keep a small mean (LRG bin 1 +0.010 σ_v, BGS −0.002 σ_v): a pre-existing ≤ 0.2σ leakage
  in the current unweighted stacks.
- B. Hadzhiyska's own figure (her samples and weights): kSZ ×0.99–1.07 (LRG), ×1.07–1.10 ± 0.1
  (BGS); f_gas both weighted ×0.93–0.98 (LRG), ×0.86–0.87 ± 0.07–0.10 (BGS). Agrees within ~1σ.

**Open question — r_V.** Weighted/unweighted ratios cancel r_V, so it is not known whether
0.65 / 0.64 apply to the weighted sample; the weights up-weight satellite-rich regions, where the
reconstruction is expected to do worse (Ondaro-Mallea et al. 2026, arXiv:2607.23339). To be asked
of B. Hadzhiyska and F. Qu: which population (fibre-assigned or complete mocks) r_V was calibrated
on, and r_V for the weighted sample.

## 3. Decisions (R. H. Liu, 2026-10-23)

| # | Decision | Words |
|---|---|---|
| — | Make completeness-weighted kSZ stacks in this repo and redo the null tests | "a plan for adding the weighting to the kSZ measurements in that repo? Include in the plan redoing the null tests since that's something I want to be certain of." |
| D1 | Weights: v1.1 `WEIGHT_COMP × WEIGHT_ZFAIL` as baseline, v2 PIP (IIP) as a check | "Yes do recommended" |
| D2 | Scope: the core stacks and the null tests first; nopairs and mass bins later | "Yes." |
| D3 | Random-position null: weight the randoms by their own `WEIGHT_COMP` | "Randoms' own weight_comp" |
| D4 | The residual mean velocity of the current unweighted stacks: report as a diagnostic, do not change the fiducial | "report as diagnostic" |
| D5 | r_V: keep 0.65 / 0.64 and flag until B. Hadzhiyska and F. Qu answer | "Yes." |

D3 note: the random catalogues carry `WEIGHT_COMP` and `WEIGHT_ZFAIL` inherited from the data
object whose redshift they copy (`TARGETID_DATA`). The data weight is `WEIGHT_COMP ×
WEIGHT_ZFAIL`; whether the randoms take `WEIGHT_ZFAIL` too (≈ 1; BGS data mean 1.001) is a detail
to confirm with R. H. Liu before running N2.

## 4. Where things are

**This repo.**
- Driver: `desi/spec_Y3/make_stacks_yaml.py` (config sections `stack`, `data`, `extra`; CLI
  overrides `--field`, `--filterType`, `--filter-cut`, `--cat-names`). `stack.zeroV` subtracts
  the mean `vR` over the loaded catalogue, **before** the catalogue mask. `stack.save: false`
  reloads the saved filtering instead of filtering the map.
- Shuffle null driver: `desi/spec_Y3/make_stacks_shuffle.py` (`shuffle_seed`, `n_seeds`,
  `bootstrap`); configs in `desi/spec_Y3/configs/null/` (`*_shuffle_*.yaml`, `*_random.yaml`).
- `thumbstack.py`: estimator lists at lines ~60–75 (`cmbHit` is None here, so
  `Est = ['ksz_uniformweight']`); `loadFiltering` ~905 (reads `{filterType}_filtmap.txt`,
  `_filtmask`, `_filtnoisestddev`, `_filtarea` from `self.pathOut`); `catalogMask` ~922;
  `computeStackedProfile` ~1163, with `ksz_uniformweight` at ~1297
  (`weights = v`, `norm = std(v)/rV/Σv²`); bootstrap `SaveCovBootstrapStackedProfile` ~1739
  (`nSamples = 10000`); velocity shuffles `SaveCovVShuffleStackedProfile` ~1834. Line numbers as
  of `669bab5`.
- Configs used by the paper: `configs/LRG_dsigma_full_zeroV.yaml`,
  `configs/LRG_dsigma_lensing_full_zeroV.yaml`, `configs/BGS_dsigma_full_cigale.yaml`,
  `configs/BGS_dsigma_lensing_full_cigale.yaml` (cat_dir `.../desi/spec_Y3/BGS_CIGALE/{cat_type}_catalogues/`).

**Traps** (from earlier sessions):
- The drivers `sys.path.append` the **main tree's** absolute path, so a worktree imports main's
  `thumbstack.py` unless `PYTHONPATH=<worktree root>` is set; log `thumbstack.__file__`.
- `loadAPRadii` matches `'DSigma'` by substring; `ftheoryGaussianProfile` has no `else` branch for
  unknown filter types.
- Environment: conda `rhliu_tSZ`. The repo's `.claude/settings.json` denies `sbatch` and
  `git push`: use `salloc`.

**Data** (all under `/pscratch/sd/r/rhliu/projects/Weak_lensing/`; read only):
- Per-object products: `ksz_measurements/ACTxDESI/spec_Y3/{,CIGALE/}{field}_no_src_with_cluster_mask_zeroV/DESIY3_{LRG_z1..z4,BGS}/`
  and `spec_Y3_DSigma_lensing/...`; randoms in `spec_Y3_random/`, shuffles in `spec_Y3_shuffle/`.
- LRG catalogues (with `TARGETID`, `WEIGHT_COMP`, `WEIGHT_ZFAIL`, `VEL_LOS_RENORM`):
  `desi/spec_Y3/LRG_catalogues_zeroV/{field}_zbin{bin}_no_src_with_cluster_mask.csv`; ThumbStack
  reads the row-aligned `catalog_*.txt` in the same folder (alignment verified through the exact
  reconstruction).
- BGS CIGALE catalogue (`RA, DEC, Z, vR, LOGMSTAR, VEL_LOS_RENORM`; no weights, no `TARGETID`):
  `desi/spec_Y3/BGS_CIGALE/BGS_catalogues/full_catalog_Y3_no_src_with_cluster_mask.csv`. Its
  weights are in the side file `desi/spec_Y3/lensing/CIGALE/lens_weights/desi_BGS_bin0.fits`
  (`RA, DEC, Z, TARGETID, WEIGHT_COMP, WEIGHT_ZFAIL`), row-aligned with it (RA/Dec/z agree to
  < 10⁻⁹), made by DESIxHSC-Lensing `data_processing/make_lens_weights.py` and checked against the
  official DA2 catalogues by `TARGETID`. The BGS cap-split stacks of Figure 25 use the non-CIGALE
  catalogue, which carries the columns.
- Random catalogues (with `WEIGHT_COMP`, `WEIGHT_ZFAIL`, `FRAC_TLOBS_TILES`):
  `desi/spec_Y3_random/{LRG,BGS}_catalogues/`.
- PIP weights: `/global/cfs/cdirs/desi/survey/catalogs/DA2/LSS/loa-v1/LSScats/v2/PIP/{LRG,BGS_BRIGHT}_{NGC,SGC}_clustering.dat.fits`
  by `TARGETID` (99.5 % of BGS and 98.9 % of LRG bin 1 stacked objects found).

## 5. Plan

**0. Setup.**
- Worktree `~/projects/repos/ThumbStack_kSZ-compweight`, new branch `ksz-compweight` from `main`.
  Check `git status` of the main tree first; never change what it has checked out.
- Run everything with `PYTHONPATH=<worktree root>`; log `thumbstack.__file__`.
- Snapshot the existing `ksz_measurements/ACTxDESI/` folders (sizes, mtimes, MD5 of
  `*_measured.txt` and `cov_*` files) before and after; nothing there may change.

**1. Estimator (`thumbstack.py`).**
- New estimator `ksz_compweight`: per-object weights from the catalogue (e.g. `Catalog.wComp`),
  masked and bootstrap-resampled with `v` and `t`; velocities re-centred on the weighted mean of
  the stacked objects; `norm = σ_w / rV / Σ w v′²` (Section 2). Under velocity shuffles the
  weights stay with the galaxies and the shuffled velocities are re-centred per realisation.
- An option to choose the estimators of a run, so a weighted run computes only `ksz_compweight`
  and never rewrites `ksz_uniformweight` files.
- A read-only `filtering_dir` option: load the filtering from the production folder, write
  everything else to a new folder; refuse `save: true` with it. (Alternative: copy the four
  per-object files; never symlink, since a write would go through to the production file.)
- Unit tests (`tests/`): weights ≡ 1 equal the unweighted estimator re-centred after the mask;
  integer weights equal repeating the objects; Σ w v′ = 0.

**2. Drivers** (`make_stacks_yaml.py`, `make_stacks_shuffle.py`).
- `stack.weights: {columns: [WEIGHT_COMP, WEIGHT_ZFAIL], file: <optional side file>}`; LRG reads
  the catalogue CSV, BGS CIGALE the side file (check RA/Dec/z row by row and the length against
  `nObj`); PIP as a variant (join by `TARGETID`).
- Defaults unchanged: a production config rerun into a test folder with `save: false` must give
  byte-identical `*_measured.txt` and covariances.

**3. Configs and outputs.** `desi/spec_Y3/configs/compweight/` for LRG bins 1–4 and BGS CIGALE,
both filters, fields full, NGC and SGC. Outputs in new folders `ksz_measurements/ACTxDESI/spec_Y3_compweight/`
and `spec_Y3_DSigma_lensing_compweight/`, file names `{filterType}_ksz_compweight_*`.

**4. Known answers before production.**
- ThumbStack's weighted profiles equal `ksz_weight_check.py`'s (extend it to print the weighted
  profile itself) to ~10⁻¹², for LRG bin 1 and BGS, both filters.
- The unit tests.

**5. Weighted stacks.** Estimator-only runs with the 10,000-sample bootstrap. Time them in a smoke
test first; one interactive allocation (CLAUDE.md: interactive QOS for ≤ 3.5 h on ≤ 4 nodes).

**6. Null tests**, weighted and unweighted side by side with the same seeds. Pass criteria, fixed
before running: PTE ≥ 0.05 for each sample and filter, and no aperture's mean null above 0.3σ.
- **N1, shuffled velocities (paper Figure 23):** 1,000 realisations per sample and filter,
  estimator only; weights stay on the galaxies. χ² and PTE of the mean null against zero with the
  shuffle covariance, and the distribution of per-realisation χ².
- **N2, random positions (Figure 22):** stacks at the random catalogues' saved per-object arrays,
  weighted by the randoms' own `WEIGHT_COMP` (D3). **Dependency:** the LRG random products are
  among the stacks the 2026-10-22 audit could not reproduce (the `Mvir` bug, fixed in `669bab5`);
  DESIxHSC-Lensing `docs/fgas_integration.md` §7–8 schedules their restack (about 3 h in two
  allocations, smoke test first). Do that restack, to a new folder, before N2.
- **N3, NGC against SGC (Figure 25):** weighted cap-split stacks; χ² consistency as in the paper.
- **N4, weight controls:** weights shuffled among galaxies with velocities fixed (the
  weighted/unweighted ratio must return to 1 within errors: shows the BGS +22 % comes from the
  weights tracking environment); leakage without and with re-centring; the weighted mean
  temperature per aperture; the residual mean velocity of the unweighted stacks (D4, diagnostic).
- **N5 (optional):** stacks on tSZ- and CIB-deprojected ILC maps, the direct test of tSZ/dust
  leakage. Needs map filtering (hours) and maps that exist; check before planning.
- **N6 (optional, DESIxHSC-Lensing):** a weighted footprint-homogeneity test (paper Figure 24:
  inside against outside HSC) with `scripts/homogeneity/` and `src/ksz_reconstruct.py`.

**7. Downstream, in DESIxHSC-Lensing** (the `lensing-v2.2` worktree until it is merged, then
`main`). Scratch copies of the plot configs reading the weighted kSZ and the v2.2 lensing make
ratio NPZs for both filters; Figures 8, 9, 14 and 16 from temporary SimulationStacker configs as
in `docs/lensing_v2.2.md` §4 (no production config changed). The unweighted pair (v2.1) is the
robustness check. r_V as D5.

**8. Docs and commits.** Code, configs and a results doc here; an update of
`docs/completeness_weights.md` in DESIxHSC-Lensing. Code review before each commit and R. H. Liu's
approval of each commit; never push.

**Out of scope for now:** nopairs and mass-bin weighted stacks (D2); recalibrating r_V (needs
mocks; D5); the adoption decision and the paper.
