"""NGC against SGC for the LRG kSZ stacks with the v2 filter (DSigma_lensing), with Liam Parker's LiFT
velocities (v_lift_R5) and the official BAO-reconstruction velocities (2026-11-17; cap stacks from
configs/velocity_dr2/LRG_dsigma_lensing_caps_*.yaml, both on the same per-object cap filtering, so only the
velocities differ; velocities zeroed in each cap and z bin).

Per LRG bin:
  1. NGC against SGC for each velocity set: chi2 of the difference with C_NGC + C_SGC, each Hartlap corrected
     (DESIxHSC-Lensing src/snr.py consistency_test, as in the paper's table). The official set must reproduce the
     paper (PTE 0.807, 0.124, 0.454, 0.548). Independent of r_V (one value for both caps).
  2. LiFT against official in each cap: amplitude ratio A = t^T W P_lift / t^T W P_off (template t = P_off, W the
     inverse bootstrap cov of P_off) and the chi2 of the difference, with the paired bootstrap (identical
     resamples in both velocity runs, seeded per sample). The noisy template biases A towards the noise slope
     tr(W C_off,lift)/9; the bias, from a Gaussian Monte Carlo with the paired bootstrap covariance and A = 1,
     is printed next to each A (larger per cap than on the full sky).
  3. Did the cap difference change: D = (NGC - SGC)_lift - (NGC - SGC)_off, covariance from the paired samples of
     each cap, summed over the two caps (different objects); also summed over the four bins.
Tests 2 and 3 with r_V (a): official / 0.65 (the paper's), LiFT R5 / Liam's Abacus altMTL values (0.553, 0.623,
0.598, 0.598; centrals; his LRG3 for our bins 3 and 4); also without r_V, where A should be near r_lift / r_off.
The catalogue velocities (LiFT against official, zeroed) correlate at r = 0.81-0.89 in each cap and bin, and the
paired bootstrap samples at 0.80-0.89, so most of the CMB x velocity noise is shared and a fluctuation would
largely persist.
The object bootstrap underestimates the errors at large apertures (up to x1.35 on the random nulls,
checks/random_diagnostic.py), so the chi2 of tests 2 and 3 are upper bounds and the errors on A lower bounds.
Profiles are in map units (uK).

Recorded run: 2026-11-17, on the cap stacks of interactive job 61499189; output in
/pscratch/sd/r/rhliu/projects/Weak_lensing/figs/vel_dr2/compare_velocity_caps.txt and compare_velocity_caps.{png,pdf},
since moved to desi/figures/2026-11/11-17/.

Usage: python compare_velocity_caps.py [--out DIR]   (login node; writes a figure and a summary)
"""
import argparse
import os
from datetime import date
import sys

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
from scipy import stats  # noqa: E402

sys.path.insert(0, '/global/u2/r/rhliu/projects/DESIxHSC-Lensing/src')
from snr import consistency_test  # noqa: E402

A_DIR = '/pscratch/sd/r/rhliu/projects/Weak_lensing/ksz_measurements/ACTxDESI/'
SUB = '{cap}_no_src_with_cluster_mask_zeroV/DESIY3_LRG_z{b}/'
EST = 'DSigma_lensing_ksz_uniformweight'
VELS = {'official': A_DIR + 'spec_Y3_DSigma_lensing_vel_dr2/official/',
        'LiFT R5': A_DIR + 'spec_Y3_DSigma_lensing_vel_dr2/v_lift_R5/'}
PROD = A_DIR + 'spec_Y3_DSigma_lensing/' # the paper's cap stacks
BINS = [1, 2, 3, 4]
CAPS = ['NGC', 'SGC']
ZLABEL = {1: '0.4 < z < 0.6', 2: '0.6 < z < 0.8', 3: '0.8 < z < 0.95', 4: '0.95 < z < 1.1'}
N_BOOT = 10000
PAPER_PTE = {1: 0.807, 2: 0.124, 3: 0.454, 4: 0.548}
RV = {'a': {'official': {b: 0.65 for b in BINS}, 'LiFT R5': {1: 0.553, 2: 0.623, 3: 0.598, 4: 0.598}},
      'raw': {'official': {b: 1. for b in BINS}, 'LiFT R5': {b: 1. for b in BINS}}}


def hartlap(n_s, n_r):
    """Hartlap factor for an inverse covariance from n_s samples of n_r bins."""
    return (n_s - n_r - 2) / (n_s - 1)


def load(base, cap, b):
    """Aperture radii [arcmin], profile, bootstrap cov and bootstrap samples of one stack."""
    d = base + SUB.format(cap=cap, b=b)
    m = np.loadtxt(d + EST + '_measured.txt')
    cov = np.loadtxt(d + 'cov_' + EST + '_bootstrap.txt')
    samples = np.loadtxt(d + EST + '_bootstrap_samples.txt') if os.path.exists(d + EST + '_bootstrap_samples.txt') else None
    return m[:, 0], m[:, 1], cov, samples


def chi2_pte(diff, cov, n_s):
    """chi2 (Hartlap) of a difference vector against zero, and its PTE."""
    c = hartlap(n_s, len(diff)) * diff @ np.linalg.solve(cov, diff)
    return c, stats.chi2.sf(c, len(diff))


def amplitude(p_off, s_off, p_lift, s_lift, n_mc=2000, seed=1):
    """Amplitude ratio A (template P_off), its paired-bootstrap error, and its template-noise bias.

    The bias is the mean A - 1 over Gaussian realisations of both profiles around P_off (A = 1), drawn
    with the joint paired-bootstrap covariance.
    """
    tw = p_off @ np.linalg.inv(np.cov(s_off, rowvar=False))
    amp, eamp = (tw @ p_lift) / (tw @ p_off), np.std((s_lift @ tw) / (s_off @ tw))
    n_r = len(p_off)
    draws = np.random.default_rng(seed).multivariate_normal(
        np.concatenate([p_off, p_off]), np.cov(np.hstack([s_off, s_lift]), rowvar=False), n_mc)
    a_mc = []
    for d in draws:
        t = d[:n_r] @ np.linalg.inv(np.cov(s_off, rowvar=False))
        a_mc.append((t @ d[n_r:]) / (t @ d[:n_r]))
    return amp, eamp, np.mean(a_mc) - 1


parser = argparse.ArgumentParser(description=__doc__.split('\n')[0])
# desi/figures/YYYY-MM/MM-DD/ of the run date (git-ignored)
FIG_DIR = os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', '..', 'figures', date.today().strftime('%Y-%m/%m-%d'))
parser.add_argument('--out', default=FIG_DIR, help='Output folder (default: desi/figures/YYYY-MM/MM-DD/).')
args = parser.parse_args()
os.makedirs(args.out, exist_ok=True)

data = {}
for vel, base in VELS.items():
    for cap in CAPS:
        for b in BINS:
            data[vel, cap, b] = load(base, cap, b)
for b in BINS:
    r0 = data['official', 'NGC', b][0]
    for key, d in data.items():
        if key[2] == b and not np.array_equal(d[0], r0):
            raise ValueError(f"{key}: aperture radii differ")
        if key[2] == b and (d[3] is None or d[3].shape != (N_BOOT, len(r0))):
            raise ValueError(f"{key}: bootstrap samples missing or not ({N_BOOT}, {len(r0)})")
    for cap in CAPS:
        s_off, s_lift = data['official', cap, b][3], data['LiFT R5', cap, b][3]
        pair_r = min(np.corrcoef(s_off[:, k], s_lift[:, k])[0, 1] for k in range(len(r0)))
        if pair_r < 0.5:
            raise ValueError(f"{cap} z{b}: bootstrap samples look unpaired (min correlation {pair_r:.2f})")

lines = ["1. NGC against SGC (chi2 with C_NGC + C_SGC, Hartlap; independent of r_V)"]
caps_test = {}
for b in BINS:
    # known answer: the official cap stacks are the paper's
    for cap in CAPS:
        prod = load(PROD, cap, b)
        if np.max(np.abs(prod[1] - data['official', cap, b][1])) > 1e-10 * np.max(np.abs(prod[1])):
            raise ValueError(f"official {cap} z{b} differs from the production cap stack")
    out = []
    for vel in VELS:
        res = consistency_test(data[vel, 'NGC', b][1], data[vel, 'NGC', b][2], data[vel, 'SGC', b][1],
                               data[vel, 'SGC', b][2], N_BOOT, N_BOOT)
        caps_test[vel, b] = res
        out.append(f"{vel}: chi2 = {res['chi2']:.1f} / {res['dof']}, PTE = {res['pte']:.3f} ({res['n_sigma']:.1f} sigma)")
    if abs(caps_test['official', b]['pte'] - PAPER_PTE[b]) > 5e-4:
        raise ValueError(f"z{b}: official NGC-SGC PTE {caps_test['official', b]['pte']:.3f} is not the paper's {PAPER_PTE[b]}")
    lines.append(f"  LRG z{b} ({ZLABEL[b]}): " + ';  '.join(out))

test3 = {}
for key in ('a', 'raw'):
    rv = RV[key]
    lines.append(f"\n2. LiFT against official in each cap, r_V {'(a)' if key == 'a' else 'none'} (paired bootstrap;"
                 " A [template bias from Monte Carlo])")
    if key == 'raw':
        lines.append("  expected A = r_lift / r_off (Abacus centrals, LiFT R5 / RecSym): "
                     + ", ".join(f"z{b} {RV['a']['LiFT R5'][b] / r:.3f}" for b, r in zip(BINS, (0.598, 0.647, 0.617, 0.617))))
    for b in BINS:
        out = []
        for cap in CAPS:
            _, p_off, _, s_off = data['official', cap, b]
            _, p_lift, _, s_lift = data['LiFT R5', cap, b]
            p_off, s_off = p_off / rv['official'][b], s_off / rv['official'][b]
            p_lift, s_lift = p_lift / rv['LiFT R5'][b], s_lift / rv['LiFT R5'][b]
            amp, eamp, bias = amplitude(p_off, s_off, p_lift, s_lift)
            c, pte = chi2_pte(p_lift - p_off, np.cov(s_lift - s_off, rowvar=False), N_BOOT)
            out.append(f"{cap} A = {amp:.3f} +- {eamp:.3f} [{bias:+.3f}], chi2 = {c:.1f}/{len(p_off)} (PTE {pte:.3f})")
        lines.append(f"  LRG z{b}: " + ';  '.join(out))
    lines.append(f"\n3. Change of the cap difference, (NGC - SGC)_LiFT - (NGC - SGC)_official, r_V {'(a)' if key == 'a' else 'none'}")
    chi2_all = 0.
    for b in BINS:
        diff, cov = 0., 0.
        for sign, cap in ((1, 'NGC'), (-1, 'SGC')):
            p = {v: data[v, cap, b][1] / rv[v][b] for v in VELS}
            s = {v: data[v, cap, b][3] / rv[v][b] for v in VELS}
            diff = diff + sign * (p['LiFT R5'] - p['official'])
            cov = cov + np.cov(s['LiFT R5'] - s['official'], rowvar=False)
        c, pte = chi2_pte(diff, cov, N_BOOT)
        test3[key, b] = pte
        chi2_all += c
        lines.append(f"  LRG z{b}: chi2 = {c:.1f} / {len(diff)}, PTE = {pte:.3f}")
    dof_all = len(BINS) * len(diff)
    lines.append(f"  all bins: chi2 = {chi2_all:.1f} / {dof_all}, PTE = {stats.chi2.sf(chi2_all, dof_all):.3f}")

summary = '\n'.join(lines)
print(summary)
with open(os.path.join(args.out, 'compare_velocity_caps.txt'), 'w') as f:
    f.write(__doc__ + '\n' + summary + '\n')

fig, axes = plt.subplots(2, 4, figsize=(16, 6.5), sharex=True, gridspec_kw={'height_ratios': [2, 1]})
rv = RV['a']
for j, b in enumerate(BINS):
    ax = axes[0, j]
    for (vel, c), dx in zip((('official', 'k'), ('LiFT R5', 'C1')), (-0.05, 0.05)):
        for cap, mfc, ddx in (('NGC', c, -0.02), ('SGC', 'none', 0.02)):
            r, p, cov, _ = data[vel, cap, b]
            ax.errorbar(r + dx + ddx, p / rv[vel][b], np.sqrt(np.diag(cov)) / rv[vel][b], fmt='o' if vel == 'official' else 's',
                        color=c, mfc=mfc, ms=5, capsize=2, label=f'{vel} {cap}')
        r = data[vel, 'NGC', b][0]
        d = (data[vel, 'NGC', b][1] - data[vel, 'SGC', b][1]) / rv[vel][b]
        e = np.sqrt(np.diag(data[vel, 'NGC', b][2] + data[vel, 'SGC', b][2])) / rv[vel][b]
        axes[1, j].errorbar(r + dx, d, e, fmt='o' if vel == 'official' else 's', color=c, capsize=2,
                            label=f"{vel}: PTE {caps_test[vel, b]['pte']:.3f}")
    ax.axhline(0, color='0.6', lw=0.8)
    ax.set_title(f"LRG z{b}, {ZLABEL[b]}")
    axes[1, j].set_title(f"change of NGC $-$ SGC (paired, test 3): PTE {test3['a', b]:.3f}", fontsize=9)
    axes[1, j].axhline(0, color='0.6', lw=0.8)
    axes[1, j].set_xlabel(r'$R$ [arcmin]')
    axes[1, j].legend(fontsize=7)
axes[0, 0].set_ylabel(r'kSZ profile, v2 filter / $r_V$ [$\mu$K]')
axes[1, 0].set_ylabel('NGC $-$ SGC\n(correlated across $R$ and\nbetween velocity sets)')
axes[0, 0].legend(fontsize=7, ncol=2)
fig.suptitle("NGC (filled) and SGC (open), v2 filter (DSigma_lensing): official / 0.65 and LiFT R5 / Abacus r_V;"
             " bootstrap errors (lower bounds at large R)", fontsize=11)
fig.tight_layout()
for ext in ('png', 'pdf'):
    fig.savefig(os.path.join(args.out, f'compare_velocity_caps.{ext}'), dpi=150)
print(f"wrote {args.out}")
