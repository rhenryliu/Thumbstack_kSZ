"""BGS_BRIGHT-21.35 kSZ stacks with the v2 filter (DSigma_lensing) and four velocity reconstructions (2026-11-19;
stacks from configs/velocity_dr2/BGS_BRIGHT-21.35_dsigma_lensing_*.yaml, all on the per-object filtering of the
official run, so only the velocities differ; velocities zeroed in each cap and z bin):
  official    DESI DR2 IFFT RecSym (recon_sm15_IFFT_recsym), the control;
  boryana_MG  Boryana Hadzhiyska's in-house MultiGrid RecSym, R = 12.5 Mpc/h, nmesh 1024;
  v_recsym    Liam Parker's RecSym (BGS share ksz_bgs_dr2_20261119; b = 1.5, f = 0.69, R = 15 Mpc/h);
  v_lift_R5   Liam Parker's LiFT, R = 5 Mpc/h smoothing.
BGS_BRIGHT-21.35 (the DR2 BAO sample of Liam's share) is a subset of the production BGS_BRIGHT-20.2, so these
stacks are not the paper's BGS stacks. Bins 1-3 are 0.1-0.2, 0.2-0.3, 0.3-0.4; bin 4 is the whole range 0.1-0.4
(it overlaps bins 1-3, so the all-bin sums use bins 1-3 only).

All four runs resample the same objects in each bootstrap realisation (seeded per sample), so the saved samples give
the covariance of any difference directly.
  Full sky, each set against the official one, and LiFT against v_recsym (the pair Liam's r values calibrate):
    the amplitude ratio A = t^T W P_x / t^T W P_ref (template t = P_ref, W the inverse bootstrap cov of P_ref) with
    its paired-bootstrap error, and the chi2 of the difference (Hartlap) and its PTE; with r_V and without.
  Caps (consistency test): 1. NGC against SGC for each set (chi2 with C_NGC + C_SGC, DESIxHSC-Lensing
    src/snr.py consistency_test; independent of r_V); 2. the change of the cap difference,
    (NGC - SGC)_x - (NGC - SGC)_official, with the paired samples of each cap summed over the caps.
r_V (Liam's Abacus altMTL values for 0.1 <= z < 0.4, centrals, one value for all bins): LiFT R5 0.560;
RecSym 0.591 for v_recsym, and for official and boryana_MG by choice (user decision 2026-11-19): on these objects
official correlates with v_recsym at r = 0.86 (NGC) / 0.95 (SGC) over 0.1-0.4 and 0.80-0.997 per cap and bin
(catalogue check of interactive job 61608703), and boryana_MG's r_V is uncalibrated.
The full-sky figure also shows, in grey, the production BGS_BRIGHT-20.2 (CIGALE) stack / 0.64 (the paper's r_V): a
different, larger sample (its full-range stack spans 0.01 < z < 0.5), so a plausibility check only.
Without r_V, A(v_lift_R5 / v_recsym) is expected near 0.560 / 0.591 = 0.948 (HalfDome: 0.665 / 0.656 = 1.014).
The object bootstrap underestimates the errors at large apertures (up to x1.35 on the random nulls,
checks/random_diagnostic.py), so the chi2 are upper bounds and the errors on A lower bounds. Profiles in uK.

Bin 4 caveat: its velocities are zeroed per cap over the whole range only, so per-cap means of up to +-113 km/s
remain in the sub-ranges (all four sets) and couple to the mean temperature; zeroing them per sub-bin moves the
bin-4 amplitude ratios (LiFT / official 1.02 -> 0.86, simple template ratio; checks/bgs_bin4_zeroing.py) into line
with bins 1-3. Read bins 1-3.
In bin 1 the official profile is near zero at small R, so the template ratio A is ill-conditioned there.

Recorded run: 2026-11-19, on the stacks of interactive jobs 61608703 (catalogues, official filtering) and 61611035
(the other three sets); output in desi/figures/2026-11/11-19/.

Usage: python compare_velocity_stacks_bgs.py [--out DIR]   (login node; writes figures and a summary)
"""
import argparse
import os
import sys
from datetime import date

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
from scipy import stats  # noqa: E402

sys.path.insert(0, '/global/u2/r/rhliu/projects/DESIxHSC-Lensing/src')
from snr import consistency_test  # noqa: E402

BASE = ('/pscratch/sd/r/rhliu/projects/Weak_lensing/ksz_measurements/ACTxDESI/spec_Y3_DSigma_lensing_vel_dr2/'
        'BGS_BRIGHT-21.35/{vel}/{field}_no_src_with_cluster_mask_zeroV/DESIY3_BGS_z{b}/')
EST = 'DSigma_lensing_ksz_uniformweight'
VELS = ['official', 'boryana_MG', 'v_recsym', 'v_lift_R5']
LABEL = {'official': 'DESI official (IFFT)', 'boryana_MG': 'Boryana MG R12.5', 'v_recsym': 'Liam RecSym',
         'v_lift_R5': 'LiFT R5'}
COLOR = {'official': 'k', 'boryana_MG': 'C0', 'v_recsym': 'C2', 'v_lift_R5': 'C1'}
BINS = [1, 2, 3, 4]
INDEP = [1, 2, 3] # bin 4 overlaps bins 1-3
ZLABEL = {1: '0.1 < z <= 0.2', 2: '0.2 < z <= 0.3', 3: '0.3 < z <= 0.4', 4: '0.1 < z <= 0.4 (all)'}
# plausibility overlay only (different sample, BGS_BRIGHT-20.2 CIGALE; its full-range stack spans 0.01 < z < 0.5)
PROD = ('/pscratch/sd/r/rhliu/projects/Weak_lensing/ksz_measurements/ACTxDESI/spec_Y3_DSigma_lensing/CIGALE/'
        'full_no_src_with_cluster_mask_zeroV/{name}/')
PROD_NAME = {1: 'DESIY3_BGS_z1', 2: 'DESIY3_BGS_z2', 3: 'DESIY3_BGS_z3', 4: 'DESIY3_BGS'}
RV_PAPER = 0.64
FIELDS = ['full', 'NGC', 'SGC']
N_BOOT = 10000
RV = {'r_V': {'official': 0.591, 'boryana_MG': 0.591, 'v_recsym': 0.591, 'v_lift_R5': 0.560},
      'no r_V': {v: 1. for v in VELS}}
PAIRS = [('boryana_MG', 'official'), ('v_recsym', 'official'), ('v_lift_R5', 'official'), ('v_lift_R5', 'v_recsym')]


def hartlap(n_s: int, n_r: int) -> float:
    """Hartlap factor for an inverse covariance from n_s samples of n_r bins."""
    return (n_s - n_r - 2) / (n_s - 1)


def load(vel: str, field: str, b: int) -> tuple:
    """Aperture radii [arcmin], profile, bootstrap cov and bootstrap samples of one stack."""
    d = BASE.format(vel=vel, field=field, b=b)
    m = np.loadtxt(d + EST + '_measured.txt')
    return m[:, 0], m[:, 1], np.loadtxt(d + 'cov_' + EST + '_bootstrap.txt'), np.loadtxt(d + EST + '_bootstrap_samples.txt')


def chi2_pte(diff: np.ndarray, cov: np.ndarray) -> tuple:
    """chi2 (Hartlap) of a difference vector against zero, and its PTE."""
    c = hartlap(N_BOOT, len(diff)) * diff @ np.linalg.solve(cov, diff)
    return c, stats.chi2.sf(c, len(diff))


def compare(p_ref: np.ndarray, s_ref: np.ndarray, p_x: np.ndarray, s_x: np.ndarray, n_mc: int = 2000,
            seed: int = 1) -> dict:
    """Amplitude ratio A of x to ref (template P_ref) and the difference chi2, with the paired bootstrap.

    The error on A is the half-width of the central 68% of its bootstrap values (the ratio is heavy-tailed when
    a resampled denominator nears 0). The template-noise bias is the mean A - 1 over Gaussian realisations of
    both profiles around P_ref (A = 1), drawn with the joint paired-bootstrap covariance.
    """
    winv = np.linalg.inv(np.cov(s_ref, rowvar=False))
    tw = p_ref @ winv
    lo, hi = np.percentile((s_x @ tw) / (s_ref @ tw), [16, 84])
    n_r = len(p_ref)
    draws = np.random.default_rng(seed).multivariate_normal(
        np.concatenate([p_ref, p_ref]), np.cov(np.hstack([s_ref, s_x]), rowvar=False), n_mc)
    bias = np.mean([(d[:n_r] @ winv @ d[n_r:]) / (d[:n_r] @ winv @ d[:n_r]) for d in draws]) - 1
    cov_diff = np.cov(s_x - s_ref, rowvar=False)
    c, pte = chi2_pte(p_x - p_ref, cov_diff)
    return {'amp': (tw @ p_x) / (tw @ p_ref), 'eamp': 0.5 * (hi - lo), 'bias': bias, 'chi2': c, 'pte': pte,
            'diff': p_x - p_ref, 'ediff': np.sqrt(np.diag(cov_diff))}


parser = argparse.ArgumentParser(description=__doc__.split('\n')[0])
# desi/figures/YYYY-MM/MM-DD/ of the run date (git-ignored)
FIG_DIR = os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', '..', 'figures', date.today().strftime('%Y-%m/%m-%d'))
parser.add_argument('--out', default=FIG_DIR, help='Output folder (default: desi/figures/YYYY-MM/MM-DD/).')
args = parser.parse_args()
os.makedirs(args.out, exist_ok=True)

data = {(v, f, b): load(v, f, b) for v in VELS for f in FIELDS for b in BINS}
pair_min = {}
for (v, f, b), d in data.items():
    r0, _, _, s0 = data['official', f, b]
    if not np.array_equal(d[0], r0) or d[3].shape != (N_BOOT, len(r0)):
        raise ValueError(f"{v} {f} z{b}: apertures or bootstrap samples differ from the official run")
    # paired resamples: the per-aperture correlation of the samples is far from 0 (an unpaired run gives ~0)
    pair_r = min(np.corrcoef(d[3][:, k], s0[:, k])[0, 1] for k in range(len(r0)))
    if v != 'official':
        pair_min[v, f, b] = pair_r
    if v != 'official' and pair_r < 0.5:
        raise ValueError(f"{v} {f} z{b}: bootstrap samples look unpaired (min correlation {pair_r:.2f})")

lines = ["Minimum per-aperture correlation of the bootstrap samples with the official run: "
         + ", ".join(f"{v} {min(r for (w, _, _), r in pair_min.items() if w == v):.2f}" for v in VELS[1:]),
         "Full sky, paired bootstrap: A(x / ref) +- 68% half-width [template-noise bias]; chi2 of x - ref (PTE)"]
full = {}
for key, rv in RV.items():
    lines.append(f"\n[{key}: " + ', '.join(f"{v} {rv[v]}" for v in VELS) + "]")
    for x, ref in PAIRS:
        lines.append(f"  {x} / {ref}:")
        chi2_all = 0.
        for b in BINS:
            _, p_ref, _, s_ref = data[ref, 'full', b]
            _, p_x, _, s_x = data[x, 'full', b]
            res = compare(p_ref / rv[ref], s_ref / rv[ref], p_x / rv[x], s_x / rv[x])
            full[key, x, ref, b] = res
            chi2_all += res['chi2'] if b in INDEP else 0.
            lines.append(f"    BGS z{b} ({ZLABEL[b]}): A = {res['amp']:.3f} +- {res['eamp']:.3f} [{res['bias']:+.3f}]; "
                         f"chi2 = {res['chi2']:.1f} / {len(p_ref)}, PTE = {res['pte']:.3f}")
        dof = len(INDEP) * len(p_ref)
        lines.append(f"    bins 1-3: chi2 = {chi2_all:.1f} / {dof}, PTE = {stats.chi2.sf(chi2_all, dof):.3f}")
lines.append("\nexpected A(v_lift_R5 / v_recsym) without r_V: Abacus 0.948, HalfDome 1.014 (Liam's README, centrals)")

lines.append("\nCaps 1. NGC against SGC (chi2 with C_NGC + C_SGC, Hartlap; independent of r_V)")
caps = {}
for b in BINS:
    out = []
    for v in VELS:
        res = consistency_test(data[v, 'NGC', b][1], data[v, 'NGC', b][2], data[v, 'SGC', b][1], data[v, 'SGC', b][2],
                               N_BOOT, N_BOOT)
        caps[v, b] = res
        out.append(f"{v}: chi2 = {res['chi2']:.1f}/{res['dof']}, PTE = {res['pte']:.3f}")
    lines.append(f"  BGS z{b}: " + ';  '.join(out))
lines.append("\nCaps 2. change of the cap difference, (NGC - SGC)_x - (NGC - SGC)_official, r_V as above")
change = {}
rv = RV['r_V']
for x in VELS[1:]:
    chi2_all = 0.
    out = []
    for b in BINS:
        diff, cov = 0., 0.
        for sign, cap in ((1, 'NGC'), (-1, 'SGC')):
            diff = diff + sign * (data[x, cap, b][1] / rv[x] - data['official', cap, b][1] / rv['official'])
            cov = cov + np.cov(data[x, cap, b][3] / rv[x] - data['official', cap, b][3] / rv['official'], rowvar=False)
        c, pte = chi2_pte(diff, cov)
        change[x, b] = pte
        chi2_all += c if b in INDEP else 0.
        out.append(f"z{b} chi2 = {c:.1f}/{len(diff)} (PTE {pte:.3f})")
    dof = len(INDEP) * len(diff)
    lines.append(f"  {x}: " + '; '.join(out) + f"; bins 1-3: chi2 = {chi2_all:.1f}/{dof}, PTE = {stats.chi2.sf(chi2_all, dof):.3f}")

summary = '\n'.join(lines)
print(summary)
with open(os.path.join(args.out, 'compare_velocity_stacks_bgs.txt'), 'w') as f:
    f.write(__doc__ + '\n' + summary + '\n')

# full sky: profiles / r_V, and each set minus the official one (paired errors)
rv = RV['r_V']
fig, axes = plt.subplots(2, 4, figsize=(16, 6.5), sharex=True, gridspec_kw={'height_ratios': [2, 1]})
offsets = dict(zip(VELS, (-0.06, -0.02, 0.02, 0.06)))
for j, b in enumerate(BINS):
    m = np.loadtxt(PROD.format(name=PROD_NAME[b]) + EST + '_measured.txt')
    axes[0, j].errorbar(m[:, 0] - 0.1, m[:, 1] / RV_PAPER, m[:, 2] / RV_PAPER, fmt='o', ms=3, capsize=2, color='0.65',
                        label=f'production BGS-20.2 / {RV_PAPER} (other sample)')
    for v in VELS:
        r, p, cov, _ = data[v, 'full', b]
        a_txt = '' if v == 'official' else f", A = {full['r_V', v, 'official', b]['amp']:.2f}"
        axes[0, j].errorbar(r + offsets[v], p / rv[v], np.sqrt(np.diag(cov)) / rv[v], fmt='o', ms=4, capsize=2,
                            color=COLOR[v], label=LABEL[v] + a_txt)
        if v != 'official':
            res = full['r_V', v, 'official', b]
            axes[1, j].errorbar(r + offsets[v], res['diff'], res['ediff'], fmt='o', ms=4, capsize=2, color=COLOR[v],
                                label=f"PTE {res['pte']:.3f}")
    axes[0, j].set_title(f"BGS-21.35 z{b}, {ZLABEL[b]}")
    for ax in axes[:, j]:
        ax.axhline(0, color='0.6', lw=0.8)
        ax.legend(fontsize=7)
    axes[1, j].set_xlabel(r'$R$ [arcmin]')
axes[0, 0].set_ylabel(r'kSZ profile, v2 filter / $r_V$ [$\mu$K]')
axes[1, 0].set_ylabel('set $-$ official\n(paired bootstrap;\ncorrelated across $R$)')
fig.suptitle("BGS_BRIGHT-21.35, v2 filter (DSigma_lensing); r_V: LiFT 0.560, others 0.591 (Boryana's uncalibrated);"
             " bootstrap errors (lower bounds at large R)", fontsize=11)
fig.tight_layout()
for ext in ('png', 'pdf'):
    fig.savefig(os.path.join(args.out, f'compare_velocity_stacks_bgs.{ext}'), dpi=150)

# caps: NGC - SGC for each set
fig, axes = plt.subplots(1, 4, figsize=(16, 3.8), sharex=True)
for j, b in enumerate(BINS):
    for v in VELS:
        r = data[v, 'NGC', b][0]
        d = (data[v, 'NGC', b][1] - data[v, 'SGC', b][1]) / rv[v]
        e = np.sqrt(np.diag(data[v, 'NGC', b][2] + data[v, 'SGC', b][2])) / rv[v]
        chg = '' if v == 'official' else f", change PTE {change[v, b]:.3f}"
        axes[j].errorbar(r + offsets[v], d, e, fmt='o', ms=4, capsize=2, color=COLOR[v],
                         label=f"{LABEL[v]}: PTE {caps[v, b]['pte']:.3f}{chg}")
    axes[j].axhline(0, color='0.6', lw=0.8)
    axes[j].set_title(f"BGS-21.35 z{b}, {ZLABEL[b]}")
    axes[j].set_xlabel(r'$R$ [arcmin]')
    axes[j].legend(fontsize=6.5)
axes[0].set_ylabel('NGC $-$ SGC, / $r_V$ [$\\mu$K]\n(correlated across $R$ and sets)')
fig.suptitle("NGC against SGC (consistency test), v2 filter (DSigma_lensing), BGS_BRIGHT-21.35", fontsize=11)
fig.tight_layout()
for ext in ('png', 'pdf'):
    fig.savefig(os.path.join(args.out, f'compare_velocity_caps_bgs.{ext}'), dpi=150)
print(f"wrote {args.out}")
