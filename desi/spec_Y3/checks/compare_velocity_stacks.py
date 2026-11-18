"""LRG kSZ stacks with the v2 filter (DSigma_lensing) using Liam Parker's LiFT velocities (v_lift_R5) against the
official BAO-reconstruction velocities (2026-11-17; stacks from configs/velocity_dr2/, both on the
production per-object filtering, so only the velocities differ).

Both stacks resample the same objects in each bootstrap realisation (seeded per sample), so the saved
bootstrap samples give the covariance of their difference directly. For each LRG bin and r_V choice:
  * the difference D = P_lift / rV_lift - P_off / rV_off, its chi2 over the apertures with the paired
    bootstrap covariance (Hartlap factor) and PTE;
  * the amplitude ratio A = t^T W P_lift / t^T W P_off (template t = P_off, W = inverse bootstrap cov of
    P_off), with its error from the same ratio in every bootstrap sample (expected 1 after r_V). The noisy
    template biases A towards the noise slope tr(W C_off,lift)/9; a Gaussian Monte Carlo puts the bias at
    under 0.01 in bins 1-3 and about 0.03 (0.3 sigma) in bin 4.
r_V choices (Liam's Abacus altMTL values, centrals; his LRG3 0.8-1.1 for our bins 3 and 4):
  (a) official / 0.65 (the paper's r_V) against LiFT R5 / (0.553, 0.623, 0.598, 0.598);
  (b) official / RecSym (0.598, 0.647, 0.617, 0.617) against the same LiFT R5;
  raw: both r_V = 1, where A estimates r_lift / r_off.
The object bootstrap underestimates the errors at large apertures (up to x1.35 on the random nulls,
checks/random_diagnostic.py), so the chi2 here are upper bounds and the errors on A lower bounds.
Profiles are in map units (uK): the DSigma_lensing filter is a weighted mean, not a sum over sr.

Recorded run: 2026-11-17, on the stacks of interactive job 61495516; output in
/pscratch/sd/r/rhliu/projects/Weak_lensing/figs/vel_dr2/compare_velocity_stacks.txt and compare_velocity_stacks_{a,b}.{png,pdf},
since moved to desi/figures/2026-11/11-17/.

Usage: python compare_velocity_stacks.py [--out DIR]   (login node; writes a figure and a summary)
"""
import argparse
import os
from datetime import date

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
from scipy import stats  # noqa: E402

BASE = ('/pscratch/sd/r/rhliu/projects/Weak_lensing/ksz_measurements/ACTxDESI/spec_Y3_DSigma_lensing_vel_dr2/'
        '{vel}/full_no_src_with_cluster_mask_zeroV/DESIY3_LRG_z{b}/')
EST = 'DSigma_lensing_ksz_uniformweight'
BINS = [1, 2, 3, 4]
ZLABEL = {1: '0.4 < z < 0.6', 2: '0.6 < z < 0.8', 3: '0.8 < z < 0.95', 4: '0.95 < z < 1.1'}
RV_LIFT = {1: 0.553, 2: 0.623, 3: 0.598, 4: 0.598}
CHOICES = {
    'a': ('official / 0.65 (paper)', {b: 0.65 for b in BINS}, RV_LIFT),
    'b': ('official / Abacus RecSym', {1: 0.598, 2: 0.647, 3: 0.617, 4: 0.617}, RV_LIFT),
    'raw': ('no r_V', {b: 1. for b in BINS}, {b: 1. for b in BINS}),
}


def load(vel, b):
    """Measured profile, its bootstrap samples, and the aperture radii [arcmin]."""
    d = BASE.format(vel=vel, b=b)
    m = np.loadtxt(d + EST + '_measured.txt')
    s = np.loadtxt(d + EST + '_bootstrap_samples.txt')
    return m[:, 0], m[:, 1], s


def compare(p_off, s_off, p_lift, s_lift, rv_off, rv_lift):
    """Difference chi2 and amplitude ratio for one bin, with the paired bootstrap."""
    p_off, s_off = p_off / rv_off, s_off / rv_off
    p_lift, s_lift = p_lift / rv_lift, s_lift / rv_lift
    n_s, n_r = s_off.shape
    hartlap = (n_s - n_r - 2) / (n_s - 1)
    diff, s_diff = p_lift - p_off, s_lift - s_off
    cov_diff = np.cov(s_diff, rowvar=False)
    chi2 = hartlap * diff @ np.linalg.solve(cov_diff, diff)
    w = np.linalg.inv(np.cov(s_off, rowvar=False)) # its normalisation cancels in A
    tw = p_off @ w
    amp = (tw @ p_lift) / (tw @ p_off)
    amp_s = (s_lift @ tw) / (s_off @ tw)
    return {'diff': diff, 'ediff': np.sqrt(np.diag(cov_diff)), 'chi2': chi2, 'dof': n_r,
            'pte': stats.chi2.sf(chi2, n_r), 'amp': amp, 'eamp': np.std(amp_s),
            'p_off': p_off, 'e_off': np.std(s_off, axis=0), 'p_lift': p_lift, 'e_lift': np.std(s_lift, axis=0)}


parser = argparse.ArgumentParser(description=__doc__.split('\n')[0])
# desi/figures/YYYY-MM/MM-DD/ of the run date (git-ignored)
FIG_DIR = os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', '..', 'figures', date.today().strftime('%Y-%m/%m-%d'))
parser.add_argument('--out', default=FIG_DIR, help='Output folder (default: desi/figures/YYYY-MM/MM-DD/).')
args = parser.parse_args()
os.makedirs(args.out, exist_ok=True)

data = {}
for b in BINS:
    r, p_off, s_off = load('official', b)
    r2, p_lift, s_lift = load('v_lift_R5', b)
    if not (np.array_equal(r, r2) and s_off.shape == s_lift.shape):
        raise ValueError(f"z{b}: apertures or bootstrap samples differ between the two runs")
    # paired resamples: the per-aperture correlation of the samples is high (0.82-0.89 when written)
    pair_r = min(np.corrcoef(s_off[:, k], s_lift[:, k])[0, 1] for k in range(len(r)))
    if pair_r < 0.5:
        raise ValueError(f"z{b}: bootstrap samples look unpaired (min correlation {pair_r:.2f})")
    data[b] = (r, p_off, s_off, p_lift, s_lift)

lines = []
results = {}
for key, (label, rv_off, rv_lift) in CHOICES.items():
    lines.append(f"\n({key}) {label}; LiFT R5 r_V = {[rv_lift[b] for b in BINS]}")
    chi2_tot = 0.
    for b in BINS:
        r, p_off, s_off, p_lift, s_lift = data[b]
        res = compare(p_off, s_off, p_lift, s_lift, rv_off[b], rv_lift[b])
        results[key, b] = res
        chi2_tot += res['chi2']
        lines.append(f"  LRG z{b} ({ZLABEL[b]}): A(LiFT/official) = {res['amp']:.3f} +- {res['eamp']:.3f};"
                     f" difference chi2 = {res['chi2']:.1f} / {res['dof']}, PTE = {res['pte']:.3f}")
    dof_tot = len(BINS) * len(r)
    lines.append(f"  all bins: chi2 = {chi2_tot:.1f} / {dof_tot}, PTE = {stats.chi2.sf(chi2_tot, dof_tot):.3f}")
lines.append("\nexpected raw A = r_lift / r_off (Abacus centrals, LiFT R5 / RecSym): "
             + ", ".join(f"z{b} {RV_LIFT[b] / CHOICES['b'][1][b]:.3f}" for b in BINS))
summary = '\n'.join(lines)
print(summary)
with open(os.path.join(args.out, 'compare_velocity_stacks.txt'), 'w') as f:
    f.write(__doc__ + summary + '\n')

for key in ('a', 'b'):
    fig, axes = plt.subplots(2, 4, figsize=(16, 6), sharex=True, gridspec_kw={'height_ratios': [2, 1]})
    for j, b in enumerate(BINS):
        r = data[b][0]
        res = results[key, b]
        ax = axes[0, j]
        ax.errorbar(r - 0.04, res['p_off'], res['e_off'], fmt='o', color='k', label='official (BAO RecSym)')
        ax.errorbar(r + 0.04, res['p_lift'], res['e_lift'], fmt='s', color='C1', label='LiFT R5')
        ax.axhline(0, color='0.6', lw=0.8)
        ax.set_title(f"LRG z{b}, {ZLABEL[b]}\nA = {res['amp']:.2f} $\\pm$ {res['eamp']:.2f}, "
                     f"$\\chi^2$ = {res['chi2']:.1f}/{res['dof']}")
        ax = axes[1, j]
        ax.errorbar(r, res['diff'], res['ediff'], fmt='o', color='C3')
        ax.axhline(0, color='0.6', lw=0.8)
        ax.set_xlabel(r'$R$ [arcmin]')
    axes[0, 0].set_ylabel(r'kSZ profile, v2 filter / $r_V$ [$\mu$K]')
    axes[1, 0].set_ylabel('LiFT $-$ official\n(paired bootstrap;\ncorrelated across $R$)')
    axes[0, 0].legend(fontsize=8)
    fig.suptitle(f"v2 filter (DSigma_lensing), velocities: ({key}) {CHOICES[key][0]}, LiFT R5 / Abacus r_V;"
                 " bootstrap errors (lower bounds at large R)")
    fig.tight_layout()
    for ext in ('png', 'pdf'):
        fig.savefig(os.path.join(args.out, f'compare_velocity_stacks_{key}.{ext}'), dpi=150)
print(f"wrote {args.out}")
