"""SNR of the kSZ stacks of the stellar-mass weighting grid (configs/mstarweight_grid/), v2 filter (DSigma_lensing):
every weighting 'ksz_mstar_{label}' present in the output folders against the unweighted stack of the same objects
('ksz_uniformweight'). Reads whichever estimators exist, so it can be run after each group of the grid.

  SNR       sqrt(chi2) against zero with the Hartlap-corrected bootstrap covariance (mstarweight_snr.detection).
  W/U       GLS amplitude of the weighted stack over the unweighted one, with the unweighted covariance, and its
            paired-bootstrap error (same resamples in every group: the seeds and the stacked objects are the same).
  sigma/U   bootstrap error of the weighted stack over the unweighted one, mean over apertures.
  Null      within-sample shuffled velocities (1,000): PTE of the mean, max |mean| / sigma; pass as in
            mstarweight_snr.py (PTE >= 0.05 and max |mean| <= 0.3 sigma).
  Caps      LRG z2 NGC against SGC (PTE), and the SGC apertures where the stack is negative.
  Regress   ksz_uniformweight and ksz_mstar_p0p67 against the earlier run's ksz_uniformweight and ksz_mstarweight
            (spec_Y3_DSigma_lensing_mstarweight/): the same objects and weights, so the profiles must be equal.

Writes the report and figures to desi/figures/YYYY-MM/MM-DD/ (the run date).
Usage: python checks/mstarweight_grid_snr.py  (login node; reads text files only)
Recorded run: 2026-11-20, after groups P, A and B (jobs 61703898, 61707251, 61710099);
desi/figures/2026-11/11-20/mstarweight_grid_snr.txt.
"""
import datetime
import glob
import os

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np

from mstarweight_snr import CAPS, FT, FULL, MEAN_MAX, PTE_MIN, ROOT as ROOT_M, cap_consistency, detection, fmt, out_dir

ROOT = ROOT_M.replace('spec_Y3_DSigma_lensing_mstarweight/', 'spec_Y3_DSigma_lensing_mstarweight_grid/')
U = 'ksz_uniformweight'
# the grid, in display order: M*^(2/3) and M*^(2/3) + M_min^(2/3) (lowest 10 %) first
ORDER = ['ksz_mstar_p0p67', 'ksz_mstar_p0p67_floor10', 'ksz_mstar_p0p33', 'ksz_mstar_p1p00', 'ksz_mstar_p0p33_floor10',
         'ksz_mstar_p0p67_floor05', 'ksz_mstar_p0p67_floor25', 'ksz_mstar_p1p00_floor10']
SAMPLES = [('full', n) for n in FULL] + [(cap, 'DESIY3_LRG_z2') for cap in CAPS]


def grid_dir(name, field):
    return out_dir(name, field).replace(ROOT_M, ROOT)


def load(d, est):
    prof = np.loadtxt(d + f'{FT}_{est}_measured.txt')
    s = dict(r=prof[:, 0], d=prof[:, 1], cov=np.loadtxt(d + f'cov_{FT}_{est}_bootstrap.txt'),
             boot=np.loadtxt(d + f'{FT}_{est}_bootstrap_samples.txt'),
             vmean=np.loadtxt(d + f'{FT}_{est}_vshufflemean.txt')[:, 1],
             vcov=np.loadtxt(d + f'cov_{FT}_{est}_vshuffle.txt'),
             vsamples=np.loadtxt(d + f'{FT}_{est}_vshuffle_samples.txt'))
    s['sig'] = np.sqrt(np.diag(s['cov']))
    s['det'] = detection(s['d'], s['cov'], n_resample=len(s['boot']))
    n = len(s['vsamples'])
    s['null'] = detection(s['vmean'], s['vcov'] / n, n_resample=n)
    s['null_max'] = np.max(np.abs(s['vmean']) / s['sig'])
    return s


def present(d):
    """Estimators with finished outputs (the shuffle samples are written last) in folder d, in display order."""
    found = {os.path.basename(p)[len(FT) + 1:-len('_measured.txt')] for p in glob.glob(d + f'{FT}_ksz_*_measured.txt')}
    found = {e for e in found if os.path.exists(d + f'{FT}_{e}_vshuffle_samples.txt')}
    unknown = sorted(found - set([U] + ORDER))
    if unknown:
        print(f"warning: {d} has estimators not in ORDER, left out: {unknown}")
    return [e for e in [U] + ORDER if e in found]


def main():
    date = datetime.date.today()
    fig_dir = os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', '..', 'figures',
                           date.strftime('%Y-%m'), date.strftime('%m-%d'))
    os.makedirs(fig_dir, exist_ok=True)
    lines = [f"# Stellar-mass weighting grid, {FT} (the v2 filter); outputs {ROOT}",
             "# w = M*^p (M* linear) or M*^p + M_min^p; objects with LOGMSTAR <= 0 dropped from every estimator"]
    st = {}
    for field, name in SAMPLES:
        d = grid_dir(name, field)
        ests = present(d)
        if U not in ests:
            lines.append(f"\n### {name}, {field}: no unweighted stack yet")
            continue
        st[(field, name)] = {e: load(d, e) for e in ests}
        s = st[(field, name)]
        u = s[U]
        lines.append(f"\n### {name}, {field}")
        lines.append("  R [arcmin]                     " + fmt(u['r'], '{:6.3f}'))
        ci = np.linalg.inv(u['cov'])
        for e in ests:
            x = s[e]
            ok = x['null']['pte'] >= PTE_MIN and x['null_max'] <= MEAN_MAX
            line = f"  {e:28s} SNR {x['det']['snr']:6.2f}"
            if e != U:
                amp = (u['d'] @ ci @ x['d']) / (u['d'] @ ci @ u['d'])
                amps = np.einsum('ij,jk,ik->i', u['boot'], ci, x['boot']) / np.einsum('ij,jk,ik->i', u['boot'], ci, u['boot'])
                line += (f"  SNR/U {x['det']['snr'] / u['det']['snr']:.3f}  W/U amplitude {amp:.3f} +- {np.std(amps):.3f}"
                         f"  sigma/U {np.mean(x['sig'] / u['sig']):.3f}")
            line += f"  null [{'pass' if ok else 'FAIL'}] PTE {x['null']['pte']:.3f}, max |mean| {x['null_max']:.3f} sigma"
            lines.append(line)
            lines.append("      stack / sigma                " + fmt(x['d'] / x['sig']))
        # regression against the earlier run: same objects, same weights
        prev = out_dir(name, field)
        for e, e_prev in ((U, U), ('ksz_mstar_p0p67', 'ksz_mstarweight')):
            if e in s:
                ref = np.loadtxt(prev + f'{FT}_{e_prev}_measured.txt')
                same = np.array_equal(ref[:, 1], s[e]['d']) and np.array_equal(
                    np.loadtxt(prev + f'cov_{FT}_{e_prev}_bootstrap.txt'), s[e]['cov'])
                dd = np.max(np.abs(ref[:, 1] - s[e]['d']) / s[e]['sig'])
                lines.append(f"  regression: {e} == earlier {e_prev} (profile and covariance): {same}; "
                             f"max |profile diff| / sigma {dd:.1e}")

    if all((cap, 'DESIY3_LRG_z2') in st for cap in CAPS):
        lines.append("\n## LRG z2, NGC against SGC")
        a, b = st[('NGC', 'DESIY3_LRG_z2')], st[('SGC', 'DESIY3_LRG_z2')]
        for e in [x for x in a if x in b]:
            chi2, pte = cap_consistency(a[e], b[e])
            neg = b[e]['r'][b[e]['d'] < 0]
            lines.append(f"  {e:28s} chi2 {chi2:5.1f} / {len(a[e]['d'])}, PTE {pte:.3f}; SGC negative at R = "
                         + (fmt(neg, '{:.3f}') if len(neg) else 'none') + f"; SGC at 6' {b[e]['d'][-1] / b[e]['sig'][-1]:+.2f} sigma")

    ests = [e for e in [U] + ORDER if any(e in s for s in st.values())]
    lines.append("\n## Summary: SNR (Hartlap-corrected bootstrap covariance); in brackets, SNR / unweighted")
    lines.append(f"  {'sample':24s} " + ' '.join(f"{e.replace('ksz_mstar_', '').replace('ksz_', ''):>17s}" for e in ests))
    for key, s in st.items():
        row = f"  {key[1] + ', ' + key[0]:24s} "
        for e in ests:
            row += f"{s[e]['det']['snr']:7.2f} ({s[e]['det']['snr'] / s[U]['det']['snr']:.3f}) " if e in s else f"{'-':>17s} "
        lines.append(row)
    text = '\n'.join(lines)
    print(text)
    with open(os.path.join(fig_dir, 'mstarweight_grid_snr.txt'), 'w') as f:
        f.write(text + '\n')

    # SNR relative to the unweighted stack, per sample, for every weighting present
    fig, ax = plt.subplots(figsize=(9, 4.2))
    labels = [e.replace('ksz_mstar_', '').replace('ksz_', '') for e in ests]
    for j, (key, s) in enumerate(st.items()):
        y = [s[e]['det']['snr'] / s[U]['det']['snr'] if e in s else np.nan for e in ests]
        ax.plot(range(len(ests)), y, 'o-', color=f'C{j}', label=f"{key[1].replace('DESIY3_', '')}, {key[0]}")
    ax.axhline(1, color='gray', lw=0.8)
    ax.set_xticks(range(len(ests)))
    ax.set_xticklabels(labels, rotation=30, ha='right')
    ax.set_ylabel('SNR / SNR unweighted')
    ax.set_title(r'kSZ SNR by stellar-mass weighting ($w = M_*^p$ or $M_*^p + M_{\min}^p$), v2 filter')
    ax.legend(fontsize=8, ncol=2)
    fig.tight_layout()
    fig.savefig(os.path.join(fig_dir, 'mstarweight_grid_snr.png'), dpi=150)

    # LRG z2 by cap: unweighted and the advisor's two weightings
    if all((f, 'DESIY3_LRG_z2') in st for f in ('full',) + CAPS):
        show = [e for e in (U, 'ksz_mstar_p0p67', 'ksz_mstar_p0p67_floor10')
                if all(e in st[(f, 'DESIY3_LRG_z2')] for f in ('full',) + CAPS)]
        fig, axes = plt.subplots(1, len(show), figsize=(5 * len(show), 4.2), sharey=True, squeeze=False)
        for ax, e in zip(axes[0], show):
            for j, (field, c) in enumerate((('full', 'k'), ('NGC', 'C0'), ('SGC', 'C3'))):
                x = st[(field, 'DESIY3_LRG_z2')][e]
                ax.errorbar(x['r'] + 0.05 * (j - 1), x['d'], x['sig'], fmt='o-', color=c, ms=4, capsize=2,
                            label=f"{field} (SNR {x['det']['snr']:.1f})")
            ax.axhline(0, color='gray', lw=0.8)
            ax.set_title(f"LRG z2, {e.replace('ksz_mstar_', '').replace('ksz_', '')}")
            ax.set_xlabel(r'$R$ [arcmin]')
            ax.legend(fontsize=9)
        axes[0][0].set_ylabel('kSZ stack, v2 filter (DSigma_lensing) [map units]')
        fig.tight_layout()
        fig.savefig(os.path.join(fig_dir, 'mstarweight_grid_lrg_z2_caps.png'), dpi=150)
    print(f"\nwrote {fig_dir}/mstarweight_grid_snr.txt and figures")


if __name__ == '__main__':
    main()
