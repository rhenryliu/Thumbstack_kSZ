# kSZ stacks weighted by stellar mass

*Written 2026-11-21 by the session that made these stacks (2026-11-19 to 2026-11-20). Commits
`7536cc4` (the M\*^(2/3) weighting) and `82d0e2b` (the floor and the grid of weightings).*

## 1. The request

From R. H. Liu's advisor, in two steps:

1. Weight the galaxies of the velocity-weighted kSZ stack by M\*^(2/3), compute the kSZ SNR for BGS and
   LRG bins 1–3 (bin 4 added) on the data, and stack LRG z2 separately in NGC and SGC, to see whether
   its SGC signal still goes negative at large apertures.
2. "Too much stochasticity at the low-mass end": weight by M\*^(2/3) + M_min^(2/3), with M_min the mean
   mass of the lowest 10 % of each sample, in the hope that this lowers the noise.

Everything here uses the v2 filter (`DSigma_lensing`) and the ACT-masked samples
(`no_src_with_cluster_mask`, `zeroV`): LRG z1–z4 and BGS with CIGALE masses, full sky, plus LRG z2 in
NGC and SGC.

## 2. Method

**Weights.** w = M\*^p, with M\* the linear stellar mass, 10^`LOGMSTAR` (not log M\*), from the
intermediate CSVs. The scale of w cancels in the estimator, so the mass units and factors of h do not
matter. With the floor, w = M\*^p + M_min^p, where M_min is the mean (linear) M\* of the lowest 10 % of
the full-sky catalogue of that sample, one value per LRG bin; NGC and SGC use the full-sky z2 value.
log10 M_min = 10.876, 10.936, 10.974, 10.942 (LRG z1–z4) and 10.001 (BGS)
(`checks/mstarweight_floor_values.py`).

**Estimator.** The one already used for the completeness weights (`ksz_compweight`, see
[completeness_weighted_ksz_results.md](completeness_weighted_ksz_results.md)):

  T̂_w(R) = σ_w / r_V · Σ w v′ T(R) / Σ w v′²,  v′ = v − Σ w v / Σ w,  σ_w² = Σ w v′² / Σ w,

over the stacked objects, after any resampling or shuffling. Code names: `ksz_mstarweight` (M\*^(2/3),
weights in `Catalog.wMstar`, config `stack.mstar_weights`) and, for the grid, `ksz_mstar_{label}`
(weights in `Catalog.kszWeights`, config `stack.mstar_weightings`, e.g. `ksz_mstar_p0p67_floor10`).

**Objects.** Objects with `LOGMSTAR <= 0` (failed fits: −99 for LRG, 0 for BGS CIGALE) are dropped from
every estimator of the run, the unweighted one included, so weighted and unweighted stacks are of the
same objects: 0.2 % (LRG z1, BGS) to 1.0 % (LRG z4) of each sample. The unweighted stacks are therefore
not byte-identical to production; §4 shows the difference.

**Errors and SNR.** 10,000 bootstrap resamples, the same for every estimator of a sample (seeded by
resample index), so ratios between weightings have paired errors. SNR = √χ² of the 9-aperture stack
against zero with the Hartlap-corrected bootstrap covariance, as `detection_summary` in DESIxHSC-Lensing
`src/snr.py`. 1,000 within-sample velocity shuffles for the null test.

**Filtering.** No map was refiltered: every run reads the per-object filtering of the production
stacks (`filtering_dir`, read only), and the NGC/SGC filtering of the completeness-weight cap run.

## 3. Results

### SNR

| Sample | Unweighted | M\*^(2/3) | M\*^(2/3) + M_min^(2/3) |
|---|---|---|---|
| LRG z1 | 9.97 | 10.48 (×1.051) | 10.46 (×1.049) |
| LRG z2 | 10.86 | 11.45 (×1.054) | 11.37 (×1.046) |
| LRG z3 | 8.59 | 8.76 (×1.019) | 8.79 (×1.023) |
| LRG z4 | 5.36 | 5.57 (×1.039) | 5.54 (×1.034) |
| BGS | 8.10 | 9.88 (×1.220) | 9.75 (×1.204) |
| LRG z2 NGC | 7.30 | 7.77 (×1.065) | 7.69 (×1.053) |
| LRG z2 SGC | 8.76 | 9.16 (×1.046) | 9.11 (×1.041) |

- **M\*^(2/3) raises the SNR in every sample**: by 2–7 % for the LRGs and 22 % for BGS, whose stellar
  masses span a wider range (N_eff/N of the weights 0.89–0.92 for LRG, 0.73 for BGS).
- **The floor lowers the noise but not the SNR.** It halves the extra bootstrap error that weighting
  costs (mean over apertures, relative to unweighted: LRG 1.02–1.03 instead of 1.04–1.07; BGS 1.14
  instead of 1.22), but gives up the same share of the gain in amplitude, so the SNR stays within 1 % of
  M\*^(2/3) alone (higher only in LRG z3).
- **The amplitude changes too**: weighted/unweighted (GLS amplitude with paired errors) 1.07–1.14 for
  M\*^(2/3) and 1.04–1.09 with the floor in the LRG bins; 1.50 ± 0.09 and 1.38 ± 0.07 for BGS. A
  weighted stack is a mass-weighted mean profile: it cannot go into f_gas unless the lensing is weighted
  the same way.

### LRG z2 in NGC and SGC

The SGC stack still goes negative with the weights. Unweighted, it is negative at R = 4.1′ to 6′
(−2.41σ at 6′) while NGC stays positive; with M\*^(2/3) it is negative at the same apertures, slightly
more so (−2.72σ at 6′), and with the floor −2.65σ. The NGC-against-SGC χ² of the difference
(C_NGC + C_SGC, Hartlap corrected) has PTE 0.101 unweighted, 0.087 with M\*^(2/3), 0.082 with the floor:
the caps differ as much as before. Up-weighting the massive galaxies does not reduce the SGC tail, which
suggests (provisionally: the weights change the mix only modestly, N_eff/N ≈ 0.9) that it does not come
mainly from the low-mass galaxies.

Note that √χ² does not see the sign: the SGC tail adds to its SNR (9.16 against 7.77 in NGC).

### The grid of weightings (extra)

Run because the best weighting is not known a priori: M\*^p with p = 1/3, 2/3, 1; M\*^(2/3) with the
floor at the lowest 5, 10 and 25 %; p = 1/3 and 1 with the 10 % floor. SNR relative to unweighted:

| Weighting | LRG z1–z4 | BGS | LRG z2 NGC / SGC |
|---|---|---|---|
| M\*^(1/3) | 1.019–1.039 | 1.144 | 1.038 / 1.036 |
| M\*^(2/3) | 1.019–1.054 | 1.220 | 1.065 / 1.046 |
| M\*^1 | 0.999–1.044 | 1.232 | 1.069 / 1.030 |
| M\*^(1/3) + floor 10 % | 1.014–1.027 | 1.106 | 1.023 / 1.023 |
| M\*^(2/3) + floor 5, 10, 25 % | 1.023–1.050 | 1.197–1.207 | 1.050–1.054 / 1.040–1.042 |
| M\*^1 + floor 10 % | 1.020–1.053 | 1.237 | 1.072 / 1.042 |

For LRG, M\*^(2/3) with or without a floor is as good as anything tried, and the floor level hardly
matters; BGS gains a little more from M\*^1. The top weightings differ by 1–2 % in SNR. No error on SNR
differences was computed, and the largest of eight correlated SNRs is biased high, so the grid does not
pick a winner. No weighting removes the SGC tail of LRG z2 (−2.56σ to −2.72σ at 6′; NGC/SGC PTE
0.080–0.104).

## 4. Checks

- **Known answer.** The stacks were rebuilt offline from the per-object filtering, with ThumbStack's
  catalogue mask (point-source cut, objects without a mass, outlier cut on the rest):
  - every weighted and unweighted ThumbStack stack matches the offline estimator to ≤ 5 × 10⁻¹³ relative
    (jobs 61616995, 61712037);
  - the earlier unweighted stacks are reproduced from all objects to ≤ 3 × 10⁻¹⁴;
  - the shift of the unweighted stacks from production, up to 0.28σ (LRG z4), comes entirely from the
    dropped objects. The outlier cut keeps the same objects, and the mean-velocity leakage changes by
    < 0.005σ. In this data the outlier cut is the same in either order of the cuts, so the check does not
    test that order.
- **Regression.** The grid's unweighted and M\*^(2/3) stacks are identical, profile and covariance, to
  those of the first run, in all seven samples.
- **Null: shuffled velocities.** Every weighted stack passes (PTE ≥ 0.12, mean ≤ 0.06σ). The unweighted
  stacks fail by PTE in LRG z1–z3, BGS and z2 NGC (mean ≤ 0.19σ). That is the known mean-velocity leakage
  of the unweighted estimator, which does not re-centre the velocities (D4 in
  completeness_weighted_ksz_results.md); the same samples failed there.
- **Inputs untouched.** No file of the production stacks or catalogues changed during the runs.
- **Tests.** `python -m pytest tests/` (30 tests), including known answers for both estimators, the
  floor and the masking of objects without a weight.

## 5. Where things are

| What | Where |
|---|---|
| Configs, M\*^(2/3) | `desi/spec_Y3/configs/mstarweight/` |
| Configs, the grid (floor included) | `desi/spec_Y3/configs/mstarweight_grid/`; three estimator groups, see the config headers |
| Stacks, M\*^(2/3) | `/pscratch/sd/r/rhliu/projects/Weak_lensing/ksz_measurements/ACTxDESI/spec_Y3_DSigma_lensing_mstarweight/` |
| Stacks, the grid | `.../ACTxDESI/spec_Y3_DSigma_lensing_mstarweight_grid/` |
| SNR, nulls, caps | `checks/mstarweight_snr.py` (first run), `checks/mstarweight_grid_snr.py` (grid) |
| Known answer | `checks/mstarweight_known_answer.py` (`--grid` for the grid) |
| Floor masses | `checks/mstarweight_floor_values.py` |
| Reports and figures | `desi/figures/2026-11/11-19/mstarweight_*` and `desi/figures/2026-11/11-20/mstarweight_grid_*` (git-ignored) |

Jobs, all interactive on `desi`: 61613574 (M\*^(2/3), 53 min), 61703898, 61707251, 61710099 (the grid in
three groups, 55–71 min on 2 nodes each), 61616995 and 61712037 (known answers). Logs:
`/pscratch/sd/r/rhliu/logs/mstarweight-*`.

## 6. Open points

- Whether any weighting is adopted is not decided. If one is, the lensing ΔΣ must be weighted the same way
  before an f_gas ratio is formed.
- The SGC tail of LRG z2 is unexplained. That the weighting leaves it unchanged is only weak evidence
  against the low-mass galaxies as its source (provisional, see §3).
- r_V is unchanged (0.65 LRG, 0.64 BGS) for the weighted stacks; SNRs and ratios do not depend on it.
