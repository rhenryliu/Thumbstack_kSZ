# Completeness-weighted kSZ stacks: results and null tests

*Written 2026-10-23 on branch `ksz-compweight` (worktree `~/projects/repos/ThumbStack_kSZ-compweight`;
merged into `main` on 2026-10-25, `1a9a289`),
from the plan in [completeness_weighted_ksz.md](completeness_weighted_ksz.md). All stacks use the
**v2 filter (`DSigma_lensing`)**; none uses the v1 filter (`DSigma`), except the checks of §2 that
reproduce production files. Every number below is copied from a file named in [§8](#8-files);
"σ" is the bootstrap error of the measured stack of the same sample, filter and estimator.
Whether the paper adopts completeness weighting is not decided here.*

---

## Summary

- **The estimator** `ksz_compweight` (ThumbStack): σ_w / r_V · Σ w v′ T / Σ w v′², with
  w = `WEIGHT_COMP × WEIGHT_ZFAIL` (DA2 `LSScats/v1.1` weights), v′ the velocity minus its weighted mean over the stacked
  objects, re-centred after every bootstrap resample and every velocity shuffle. With all new
  options off, the production stacks come out byte-identical; the weighted stacks equal an
  independent offline computation to ≤ 2.4 × 10⁻¹⁴.
- **Weighted / unweighted** (GLS amplitude, unweighted covariance, paired bootstrap):
  LRG z1 0.991 ± 0.038, z2 1.028 ± 0.037, z3 1.053 ± 0.046, z4 1.047 ± 0.067, BGS 1.219 ± 0.077.
  PIP (IIP) weights agree: 0.983, 1.015, 1.079, 0.963, 1.214.
- **Weight control:** with the weights shuffled among the galaxies the ratio returns to 1.00 for
  every sample; the BGS +22 % is 2.9σ from that, so it comes from the weights tracking environment.
- **Null tests of the weighted stacks** (pass: PTE ≥ 0.05 and, for shuffle means, no aperture
  above 0.3σ): 18 of 20 pass (10 in N1, 5 in N2, 5 in N3; §4). The two failures are marginal
  shuffle nulls (PTE 0.038 and 0.037, largest aperture mean 0.07σ and 0.09σ), and an independent
  repeat of each passes (PTE 0.71 and 0.54). No aperture mean exceeds 0.09σ in any shuffle test.
- **Random positions:** the weighted nulls pass in every sample, with the randoms weighted by their
  own `WEIGHT_COMP × WEIGHT_ZFAIL` or `WEIGHT_COMP` alone (the two agree to ≤ 0.02σ). The unweighted
  LRG z2 random null fails (PTE 0.007, every aperture −1.2σ to −3.2σ). Not the mean-velocity
  leakage (≤ 0.07σ for randoms); with a 200-patch sky jackknife, whose errors are ×1.0–1.35 the
  object bootstrap's, its PTE is 0.052. **The object bootstrap underestimates the random-null errors
  at large apertures; whether it does so for the data stacks too was not checked.**
- **D4, the unweighted stacks:** the velocities of the stacked objects keep a small mean
  (+0.010 σ_v for LRG z1, −0.002 σ_v for BGS), which leaks the mean temperature into the stack by up
  to +0.19σ (LRG z1, 1′). Shuffling velocities within the stacked objects shows exactly this
  offset, so that unweighted null fails for LRG z1–z3 and BGS; the whole-catalogue shuffle of the
  published paper's Figure 23 does not see it. The weighted estimator removes it by construction.
  Reported as a diagnostic; the unweighted fiducial is unchanged (your decision D4).
- **r_V** is 0.65 / 0.64 for the weighted stacks too, unchanged and flagged (D5): every ratio here
  cancels it.

## 1. What was run

| Item | Config (`desi/spec_Y3/configs/compweight/`) | Output (`ksz_measurements/ACTxDESI/`) |
|---|---|---|
| Weighted and unweighted stacks, full sky, 10,000 bootstraps, 1,000 within-sample shuffles | `LRG_dsigma_lensing_full.yaml`, `BGS_dsigma_lensing_full_cigale.yaml` | `spec_Y3_DSigma_lensing_compweight/{,CIGALE/}full_no_src_with_cluster_mask_zeroV/` |
| Same with `LSScats/v2/PIP` (IIP) weights | `*_full*_pip.yaml` | `spec_Y3_DSigma_lensing_compweight_pip/` |
| Whole-catalogue shuffle null, 1,000 seeds (master seed 100; the first 100 are the seeds of the published Figure 23) | `*_dsigma_lensing_shuffle*.yaml` | `spec_Y3_DSigma_lensing_compweight_shuffle/` |
| NGC and SGC, new v2-filter filtering | `*_dsigma_lensing_caps.yaml --field NGC/SGC` | `spec_Y3_DSigma_lensing_compweight/{NGC,SGC}_no_src_with_cluster_mask_zeroV/` |
| Random positions, new v2-filter filtering, randoms' `WEIGHT_COMP × WEIGHT_ZFAIL` | `*_dsigma_lensing_random.yaml` | `spec_Y3_DSigma_lensing_compweight_random/` |
| Same, randoms' `WEIGHT_COMP` only | `*_dsigma_lensing_random_comp.yaml` | `spec_Y3_DSigma_lensing_compweight_random_comp/` |

Samples: LRG z1–z4 (ACT-masked, `LRG_catalogues_zeroV`), BGS CIGALE (full sky; weights from the
lensing side file `lensing/CIGALE/lens_weights/desi_BGS_bin0.fits`). The caps use the non-CIGALE BGS
catalogue, as Figure 25 did (no CIGALE cap catalogues exist). Unweighted and weighted stacks of one
run use the same objects, bootstrap resamples and shuffles, so their ratios have paired errors.

Jobs (all interactive, `desi`, except 60177318): 60171703 (checks, 1h29), 60176230 (random smoke
test, shuffles, LRG stacks, 1h11), 60178685 (PIP, LRG z1 randoms, N4), 60180989 (caps),
60182679 (LRG z2, z3 randoms), 60185347 (random `WEIGHT_COMP` variant, BGS cap check),
60187823 and 60188069 (random diagnostics); 60177318 (regular, BGS randoms, submitted by you;
5h04, peak memory 397 GB of the node's 476 GB, so a larger random catalogue would not fit).

## 2. Checks before the results

- **Defaults unchanged.** The production configs rerun with the new code, reading the production
  filtering, give `*_ksz_uniformweight_measured.txt` and `cov_*_bootstrap.txt` byte-identical to
  production for LRG z1 and BGS CIGALE with both filters, and the theory files too, except
  `DSigma_ksz_uniformweight_theory_tsz.txt` of LRG z1: production's has ~10⁻³¹² values in its profile
  column, the uninitialised `integratedY` of the drivers before `579f022`/`669bab5`; the rerun has
  zeros. Every unweighted v2-filter stack of the new folders (LRG z1–z4, BGS) is byte-identical
  to production.
- **Known answer.** `ksz_compweight` equals the weighted estimator of
  `DESIxHSC-Lensing-v2.2/scratch/v2.2_ksz_weights/ksz_weight_check.py`, recomputed offline from the
  production per-object files, to 5 × 10⁻¹⁵ – 2.4 × 10⁻¹⁴ (LRG z1–z4, BGS). The LRG z1 and BGS ratios
  and PIP ratios equal that script's tables.
- **Shuffle driver.** `make_stacks_shuffle.py` with `load_once` reproduces the published Figure 23
  shuffle files byte for byte (LRG z1–z4, BGS; 100 seeds, master seed 100).
  **Stale config found:** `configs/null/LRG_dsigma_shuffle_zeroV.yaml` says master seed 200, but the
  files were made with 100 (seed 200 does not reproduce them; the 2026-10-22 Mvir audit found the
  same). Not changed; your call.
- **Cap filtering.** Every NGC and SGC object's new per-object filtering (all four files) equals its
  row in the production full-sky v2-filter files exactly (LRG z1–z4; BGS: the 98.2 % (NGC) and
  99.5 % (SGC) of the non-CIGALE cap objects that are in the CIGALE full-sky catalogue).
- **Unit tests:** `tests/test_ksz_compweight.py`, 13 tests (22 with `test_dsigma_lensing.py`), pass.

## 3. Weighted and unweighted stacks

GLS amplitude of weighted over unweighted, with the unweighted bootstrap covariance; error from the
10,000 paired bootstrap resamples. Per aperture (1′ to 6′) in `nulls_summary.txt`.

| Sample | W/U, `LSScats/v1.1` weights | W/U, PIP weights (objects in PIP) | `LSScats/v1.1` weight over the catalogue: mean, N_eff/N |
|---|---|---|---|
| LRG z1 | 0.991 ± 0.038 | 0.983 ± 0.052 (99.0 %) | 1.099, 0.895 |
| LRG z2 | 1.028 ± 0.037 | 1.015 ± 0.066 (99.0 %) | 1.107, 0.887 |
| LRG z3 | 1.053 ± 0.046 | 1.079 ± 0.067 (99.0 %) | 1.111, 0.887 |
| LRG z4 | 1.047 ± 0.067 | 0.963 ± 0.121 (99.0 %) | 1.097, 0.905 |
| BGS | 1.219 ± 0.077 | 1.214 ± 0.074 (99.6 %) | 1.273, 0.757 |

- The weighted stacks are noisier: their bootstrap errors are ×1.06–1.10 (LRG) and ×1.15–1.17 (BGS)
  the unweighted ones.
- The BGS ratio is flat in aperture (1.22–1.38 per aperture, errors 0.10–0.21).

## 4. Null tests

Criteria fixed before the runs: PTE ≥ 0.05 for each sample and estimator, and for the means over
shuffled-velocity realisations no aperture above 0.3σ. (I applied the 0.3σ criterion to the shuffle
means only: a single random-position stack's own error is 0.39–0.44σ per aperture.)

### N1, shuffled velocities (1,000 realisations each)

(a) whole catalogue shuffled before the catalogue mask, as for the published Figure 23;
(b) velocities shuffled among the stacked objects (ThumbStack's own shuffle). Weights stay on the
galaxies; the weighted estimator re-centres each realisation.

| Sample | Weighted (a) | Weighted (b) | Unweighted (a) | Unweighted (b) |
|---|---|---|---|---|
| LRG z1 | 0.298 | **0.038** | 0.247 | **0.000** |
| LRG z2 | 0.835 | 0.663 | 0.732 | **0.000** |
| LRG z3 | **0.037** | 0.373 | **0.022** | **0.000** |
| LRG z4 | 0.777 | 0.371 | 0.843 | 0.321 |
| BGS | 0.335 | 0.489 | 0.130 | **0.000** |

PTE of the mean null against zero (bold: below 0.05). Largest aperture mean: ≤ 0.086σ weighted,
≤ 0.173σ unweighted (b). The per-realisation χ² follow χ²(9) in every case (mean 9.00–9.03; KS
p 0.08–0.96).

- **The two weighted failures**, by the fixed criteria, are LRG z1 (b) and LRG z3 (a). Repeats
  with independent realisations, run as diagnostics after seeing them (they do not replace the
  results above): LRG z1 (b), 2,000 new within-sample shuffles, PTE 0.712; LRG z3 (a), 1,000 new
  seeds (master seed 300), weighted PTE 0.535, unweighted 0.719. The LRG z3 (a) pattern is shared by
  the unweighted null with the same seeds. Both look like fluctuations: 2 failures in 10 weighted
  tests has a probability of about 9 % by chance.
- **The unweighted (b) failures are D4**, quantified in N4: the mean of the unweighted (b) null
  equals, aperture by aperture, the predicted leakage plus the noise shared with the weighted null
  (same permutations).

### N2, random positions

Random catalogues `desi/spec_Y3_random/{LRG,BGS}_catalogues/` (LRG z1–z4: 3.18M, 4.86M, 3.88M,
1.58M; BGS 14.1M, non-CIGALE), new v2-filter filtering; the randoms carry the velocity field and
the weights of the data object whose redshift they copy. χ² of each stack against zero with its own
10,000-sample bootstrap covariance (Hartlap):

| Sample | Weighted, `WEIGHT_COMP × WEIGHT_ZFAIL` | Weighted, `WEIGHT_COMP` | Unweighted |
|---|---|---|---|
| LRG z1 | 0.684 | 0.683 | 0.864 |
| LRG z2 | 0.155 | 0.157 | **0.007** |
| LRG z3 | 0.866 | 0.865 | 0.608 |
| LRG z4 | 0.180 | 0.182 | 0.131 |
| BGS | 0.247 | 0.249 | 0.090 |

- **`WEIGHT_ZFAIL` on the randoms** (your D3 detail) makes no difference: the two weightings agree
  to ≤ 0.02σ per aperture.
- **The unweighted LRG z2 failure**, diagnosed afterwards (`random_diagnostic.txt`; the result above
  stands as recorded): every aperture is negative, −1.2σ to −3.2σ, and the weighted stack has the
  same shape at −1.1 to −2.5 times the unweighted stack's error, so it is not specific to the weights.
  - Not the mean-velocity leakage: the stacked randoms' velocities have mean +0.013 σ_v and the mean
    filtered temperature at random positions is −0.006 to −0.032 map units, so the leakage is
    −0.03σ to −0.07σ for LRG z2 and ≤ 0.06σ in the other samples.
  - The errors: randoms outnumber the galaxies 5–6 times and oversample the map, so neighbouring
    apertures share pixels and the object bootstrap treats correlated objects as independent. A
    delete-one jackknife over 200 contiguous sky patches gives errors ×1.0–1.35 the bootstrap's,
    growing with aperture (LRG z2 ×1.16 at 6′, LRG z3 ×1.35). With the jackknife covariance the
    PTEs are: LRG z2 unweighted 0.052, re-centred 0.058, weighted 0.335; all others ≥ 0.08.
  - The data stacks have fewer objects per area, but their apertures overlap too; their bootstrap
    errors were not compared with a jackknife here.

### N3, NGC against SGC

χ² of the difference with C_NGC + C_SGC, each Hartlap corrected, as `src/snr.py consistency_test`
in DESIxHSC-Lensing. All pass:

| Sample | Weighted | Unweighted |
|---|---|---|
| LRG z1 | 0.594 | 0.807 |
| LRG z2 | 0.089 | 0.124 |
| LRG z3 | 0.648 | 0.454 |
| LRG z4 | 0.657 | 0.548 |
| BGS (non-CIGALE) | 0.292 | 0.165 |

### N4, weight controls (`n4_weight_controls.txt`)

| Sample | W/U, weights shuffled among galaxies (200) | Observed W/U from it | Leakage without re-centring, largest | D4: mean velocity of the stacked objects | D4: leakage in the unweighted stack, largest |
|---|---|---|---|---|---|
| LRG z1 | 1.002 ± 0.036 | −0.3σ | +0.22σ (1′) | +0.0104 σ_v | +0.19σ (1′) |
| LRG z2 | 1.001 ± 0.033 | +0.8σ | +0.05σ (1′) | +0.0048 σ_v | +0.08σ (1′) |
| LRG z3 | 1.000 ± 0.044 | +1.2σ | −0.19σ (2.9′) | +0.0112 σ_v | +0.12σ (1′) |
| LRG z4 | 0.997 ± 0.054 | +0.9σ | +0.02σ | −0.0006 σ_v | +0.005σ |
| BGS | 1.009 ± 0.072 | +2.9σ | −0.17σ (1′) | −0.0023 σ_v | −0.11σ (1′) |

- "Leakage without re-centring": the weighted stack with the velocities not re-centred on their
  weighted mean, minus the weighted stack. Re-centring removes it exactly.
- The weighted mean temperature is lower than the unweighted one at every aperture beyond 2′ (more
  tSZ in crowded regions, e.g. BGS at 3.5′: −0.027 → −0.048 map units); it enters only through the
  mean velocity, which the weighted estimator sets to zero.
- D4 leakage = ThumbStack's unweighted stack minus the same stack with the velocities re-centred
  after the mask; it equals (⟨v⟩/σ_v) ⟨T⟩ / σ to the third decimal.

## 5. Deviations from the plan, and choices made

- **v2 filter only**, your instruction; the v1 filter (`DSigma`) ran only to reproduce production
  files (§2). So the LRG random restack of the plan was made with the v2 filter (new filtering of
  all 13.5M LRG and 14.1M BGS randoms), and the v1-filter LRG random products remain unreproduced.
- **Two shuffle nulls**, 1,000 realisations each: the whole-catalogue shuffle of the published
  Figure 23 (your request; master seed 100, so its first 100 seeds are Figure 23's) and the
  within-sample shuffle, which sees the D4 leakage. 1,000 seeds were affordable through the new
  `load_once` option.
- **Bootstrap:** the weighted estimator re-centres the velocities on each resample (a bootstrap of
  the whole estimator); `ksz_weight_check.py` used the full-sample mean. The two differ only through
  the resampled mean velocity, ~σ_v/√N; not measured.
- **PIP:** objects without a PIP weight (0.4–1.0 %) are dropped from both estimators of the PIP run.
- **Randoms' weights:** both `WEIGHT_COMP × WEIGHT_ZFAIL` (like the data) and `WEIGHT_COMP` alone,
  your choice; they agree.
- **The 0.3σ criterion** was applied to the means over shuffle realisations only: a single random
  stack's own bootstrap error is 0.39–0.44σ per aperture (weighted, all samples), so the criterion
  would fail on noise alone.
- **Diagnostics run after seeing failures** (independent shuffle repeats, the random jackknife) are
  labelled as such and do not replace the pre-set results.
- **Jobs:** interactive allocations of under 2 h (the limit of the session's shell) and one regular
  job for the BGS randoms, which you submitted.

## 6. Not done, open

- **N5** (deprojected ILC maps) and **N6** (footprint homogeneity): not run.
- **Bootstrap against jackknife for the data stacks** (N2 diagnostic): not checked; a separate
  decision, since it would change the errors of the published unweighted stacks too.
- **v1 filter (`DSigma`):** configs exist (`*_dsigma_*.yaml`, outputs in `spec_Y3_DSigma_compweight*`)
  but were not run, per your instruction; run them only if a v1-filter result is asked for.
- **Downstream f_gas** with lensing v2.2 (plan step 7, DESIxHSC-Lensing): done after this doc was
  first committed, on 2026-10-24, as one-off copies, with the v2 filter only. f_gas v2.2 (both
  probes weighted) against v2.1 (neither) in Figures 8, 9, 14 and 16: SimulationStacker's git-ignored
  `figures/2026-10/10-24/completeness_weights_DSigma_lensing/` (with `completeness_weights_summary.txt`
  and `fgas_measurements_v2.1_v2.2.csv`); ratio NPZs in `ratio_measurements/spec_Y3_v2.2_compweight/`
  and `ratio_measurements/compweight_v2.1_v2.2/`; scripts in DESIxHSC-Lensing's git-ignored
  `scratch/compweight_figs/`. The paper-figure runner does not read them.
- **r_V for the weighted sample** (D5): open, for B. Hadzhiyska and F. Qu.
- Nopairs and mass-bin weighted stacks: out of scope (D2); `apply_catalog_mask` does not carry the
  weights, so the drivers refuse weights together with a mass filter.

## 7. Code

`f013ac4` (ThumbStack_kSZ, branch `ksz-compweight`): `ksz_compweight` and the options `estimators`,
`filteringDir` (read-only per-object filtering), `nVShuffle`, `saveSamples` in `thumbstack.py`;
`stack.weights`, `estimators`, `filtering_dir`, `vshuffle`, `n_vshuffle`, `save_samples` in the
drivers, `load_once` in `make_stacks_shuffle.py`; `desi/spec_Y3/weight_utils.py`; tests. Summary of
the null tests: `desi/spec_Y3/compweight_nulls.py`. One-off checks, committed 2026-10-25 from the
worktree's `scratch/`: `desi/spec_Y3/checks/` (`compweight_known_answer.py`, `compweight_n4.py`,
`cap_filtering_check.py`, `random_diagnostic.py`) and the LRG z3 seed-300 repeat config
`desi/spec_Y3/configs/compweight/checks/shuffle_seed300_LRG_dsigma_lensing.yaml`.

## 8. Files

Under `/pscratch/sd/r/rhliu/projects/Weak_lensing/ksz_measurements/ACTxDESI/`:

| File | Contents |
|---|---|
| `spec_Y3_DSigma_lensing_compweight/nulls_summary.txt` | `compweight_nulls.py` output: §3, N1–N3 |
| `spec_Y3_DSigma_lensing_compweight/n4_weight_controls.txt` | N4, the independent within-sample shuffles |
| `spec_Y3_DSigma_lensing_compweight/known_answer_LRG_z2-z4.txt` | §2 known answer, LRG z2–z4 (z1 and BGS: log of job 60171703) |
| `spec_Y3_DSigma_lensing_compweight/cap_filtering_check_LRG.txt` | §2 cap filtering |
| `compweight_checks/shuffle_seed300/summary.txt` | N1 repeat for LRG z3 (a) |
| `spec_Y3_DSigma_lensing_compweight_random/random_diagnostic.txt` | N2 leakage and jackknife |
| `compweight_checks/defaults/`, `compweight_checks/shuffle_known_answer_DSigma/` | §2 checks against production |

Logs: `/pscratch/sd/r/rhliu/logs/compweight-A*-*.log`, `cw-bgs-random-60177318.{out,err}`.
