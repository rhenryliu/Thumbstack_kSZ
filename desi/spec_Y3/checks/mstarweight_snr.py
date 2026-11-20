"""SNR of the stellar-mass-weighted kSZ stacks (w = M*^(2/3), estimator 'ksz_mstarweight') against the
unweighted stacks of the same objects ('ksz_uniformweight'), v2 filter (DSigma_lensing), from the outputs of
configs/mstarweight/ (one-off check, kept as the record of that run).

  SNR       sqrt(chi2) of the stack against zero with the Hartlap-corrected bootstrap covariance (its
            number of resamples, 10,000), PTE and the Gaussian-equivalent n_sigma, as DESIxHSC-Lensing src/snr.py
            detection_summary. LRG z1-z4 and BGS (CIGALE) full sky; LRG z2 NGC and SGC.
  M/U       weighted / unweighted per aperture, with the linearised paired-bootstrap error (same
            resamples), and as a GLS amplitude with the unweighted covariance.
  Known     the unweighted stack here against the earlier one of the same sample (production full sky,
            completeness-weight run for the caps): it differs only by the objects with LOGMSTAR <= 0.
  Null      within-sample shuffled velocities (1,000): chi2 of the mean against zero with C/N (Hartlap),
            and max |mean| / bootstrap sigma. Pass: PTE >= 0.05 and max |mean| <= 0.3 sigma, as for the
            completeness-weighted stacks (compweight_nulls.py).
  Caps      LRG z2, NGC against SGC: chi2 of the difference with C_NGC + C_SGC, each Hartlap corrected
            (src/snr.py consistency_test).

Writes the report and two figures to desi/figures/YYYY-MM/MM-DD/ (the run date).
Usage: python checks/mstarweight_snr.py  (login node; reads text files only)
"""
import datetime
import os
import sys

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
from scipy import special, stats

K = '/pscratch/sd/r/rhliu/projects/Weak_lensing/ksz_measurements/ACTxDESI/'
ROOT = K + 'spec_Y3_DSigma_lensing_mstarweight/'
FT = 'DSigma_lensing'
CUT = 'no_src_with_cluster_mask_zeroV'
EST = {'U': 'ksz_uniformweight', 'M': 'ksz_mstarweight'}
PTE_MIN, MEAN_MAX = 0.05, 0.3
FULL = [f'DESIY3_LRG_z{b}' for b in (1, 2, 3, 4)] + ['DESIY3_BGS']
CAPS = ('NGC', 'SGC')


def out_dir(name, field='full'):
    """Folder of one stack in this run; the BGS stacks use the CIGALE catalogue."""
    sub = 'CIGALE/' if name == 'DESIY3_BGS' else ''
    return f'{ROOT}{sub}{field}_{CUT}/{name}/'


def earlier_dir(name, field='full'):
    """The earlier unweighted stack of the same sample (the known answer)."""
    if field == 'full':
        sub = 'CIGALE/' if name == 'DESIY3_BGS' else ''
        return f'{K}spec_Y3_DSigma_lensing/{sub}full_{CUT}/{name}/'
    return f'{K}spec_Y3_DSigma_lensing_compweight/{field}_{CUT}/{name}/'


def hartlap(n, p):
    return (n - p - 2) / (n - 1)


def detection(d, cov, n_resample):
    """sqrt(chi2) against zero with the Hartlap-corrected covariance, PTE and n_sigma."""
    chi2 = float(d @ np.linalg.inv(cov / hartlap(n_resample, len(d))) @ d)
    pte = float(stats.chi2.sf(chi2, df=len(d)))
    return dict(snr=np.sqrt(chi2), chi2=chi2, pte=pte, n_sigma=float(np.sqrt(2) * special.erfcinv(pte)))


def load_stack(d):
    """Measured profile, bootstrap covariance and samples, and shuffle products of both estimators."""
    out = {}
    for k, e in EST.items():
        prof = np.loadtxt(d + f'{FT}_{e}_measured.txt')
        out[k] = dict(r=prof[:, 0], d=prof[:, 1],
                      cov=np.loadtxt(d + f'cov_{FT}_{e}_bootstrap.txt'),
                      boot=np.loadtxt(d + f'{FT}_{e}_bootstrap_samples.txt'),
                      vmean=np.loadtxt(d + f'{FT}_{e}_vshufflemean.txt')[:, 1],
                      vcov=np.loadtxt(d + f'cov_{FT}_{e}_vshuffle.txt'),
                      vsamples=np.loadtxt(d + f'{FT}_{e}_vshuffle_samples.txt'))
        out[k]['sig'] = np.sqrt(np.diag(out[k]['cov']))
    return out


def fmt(x, f='{:+.3f}'):
    return ' '.join(f.format(v) for v in x)


def report_sample(label, s, earlier, lines):
    u, m = s['U'], s['M']
    if np.array_equal(u['d'], m['d']):
        raise ValueError(f"{label}: weighted and unweighted stacks are identical; the weights did not act")
    lines.append(f"\n### {label}")
    lines.append("  R [arcmin]          " + fmt(u['r'], '{:6.3f}'))
    for k in ('U', 'M'):
        lines.append(f"  {k} stack / sigma     " + fmt(s[k]['d'] / s[k]['sig']))
    for k in ('U', 'M'):
        det = detection(s[k]['d'], s[k]['cov'], n_resample=len(s[k]['boot']))
        s[k]['det'] = det
        lines.append(f"  {k} {EST[k]}: SNR {det['snr']:.2f} (chi2 {det['chi2']:.1f} / {len(u['d'])}, "
                     f"PTE {det['pte']:.2e}, n_sigma {det['n_sigma']:.2f})")
    ratio = m['d'] / u['d']
    # linearised (delta-method) paired error: the ratio of resamples is heavy tailed where U is noisy
    err = np.std((m['boot'] - ratio * u['boot']) / u['d'], axis=0)
    ci = np.linalg.inv(u['cov'])
    amp = (u['d'] @ ci @ m['d']) / (u['d'] @ ci @ u['d'])
    amps = np.einsum('ij,jk,ik->i', u['boot'], ci, m['boot']) / np.einsum('ij,jk,ik->i', u['boot'], ci, u['boot'])
    lines.append("  M/U per aperture    " + fmt(ratio, '{:6.3f}'))
    lines.append("    paired error      " + fmt(err, '{:6.3f}'))
    lines.append(f"  M/U amplitude (GLS, unweighted covariance) {amp:.3f} +- {np.std(amps):.3f}")
    lines.append("  bootstrap sigma M/U " + fmt(m['sig'] / u['sig'], '{:6.3f}'))
    # known answer: the earlier unweighted stack, on 0.2-1 % more objects
    ref = np.loadtxt(earlier + f'{FT}_ksz_uniformweight_measured.txt')
    if not np.allclose(ref[:, 0], u['r']):
        raise ValueError(f"{label}: aperture radii differ from {earlier}")
    lines.append("  known answer, (U - earlier U) / sigma_U " + fmt((u['d'] - ref[:, 1]) / u['sig']))
    # null: velocities shuffled among the stacked objects
    for k in ('U', 'M'):
        n = len(s[k]['vsamples'])
        null = detection(s[k]['vmean'], s[k]['vcov'] / n, n_resample=n)
        worst = np.max(np.abs(s[k]['vmean']) / s[k]['sig'])
        ok = null['pte'] >= PTE_MIN and worst <= MEAN_MAX
        lines.append(f"  null, {k} shuffled velocities ({n}): [{'pass' if ok else 'FAIL'}] chi2 {null['chi2']:.1f} / "
                     f"{len(u['d'])}, PTE {null['pte']:.3f}; max |mean| {worst:.3f} sigma")
        lines.append("      mean / sigma " + fmt(s[k]['vmean'] / s[k]['sig']))
    return ratio, err


def cap_consistency(a, b):
    """chi2 of NGC - SGC with C_NGC + C_SGC, each Hartlap corrected; PTE."""
    delta = a['d'] - b['d']
    cov = a['cov'] / hartlap(len(a['boot']), len(delta)) + b['cov'] / hartlap(len(b['boot']), len(delta))
    chi2 = float(delta @ np.linalg.inv(cov) @ delta)
    return chi2, float(stats.chi2.sf(chi2, df=len(delta)))


def main():
    date = datetime.date.today()
    fig_dir = os.path.join(os.path.dirname(__file__), '..', '..', 'figures', date.strftime('%Y-%m'), date.strftime('%m-%d'))
    os.makedirs(fig_dir, exist_ok=True)
    lines = [f"# kSZ weighted by M*^(2/3) (ksz_mstarweight, M) vs unweighted (ksz_uniformweight, U), {FT} (the v2 filter)",
             f"# outputs: {ROOT}", "# objects with LOGMSTAR <= 0 dropped from both M and U"]
    stacks = {}
    for name in FULL:
        stacks[('full', name)] = load_stack(out_dir(name))
        report_sample(f"{name}, full sky", stacks[('full', name)], earlier_dir(name), lines)
    for cap in CAPS:
        stacks[(cap, 'DESIY3_LRG_z2')] = load_stack(out_dir('DESIY3_LRG_z2', cap))
        report_sample(f"DESIY3_LRG_z2, {cap}", stacks[(cap, 'DESIY3_LRG_z2')],
                      earlier_dir('DESIY3_LRG_z2', cap), lines)

    lines.append("\n## LRG z2, NGC against SGC")
    for k in ('U', 'M'):
        chi2, pte = cap_consistency(stacks[('NGC', 'DESIY3_LRG_z2')][k], stacks[('SGC', 'DESIY3_LRG_z2')][k])
        dof = len(stacks[('SGC', 'DESIY3_LRG_z2')][k]['d'])
        lines.append(f"  {k} {EST[k]}: chi2 {chi2:.1f} / {dof}, PTE {pte:.3f}")
    sgc = stacks[('SGC', 'DESIY3_LRG_z2')]
    for k in ('U', 'M'):
        neg = sgc[k]['r'][sgc[k]['d'] < 0]
        lines.append(f"  SGC {k}: negative at R = " + (fmt(neg, '{:.3f}') if len(neg) else 'none') + " arcmin")

    lines.append("\n## Summary: SNR (Hartlap-corrected bootstrap covariance)")
    lines.append(f"  {'sample':28s} {'SNR U':>7s} {'SNR M':>7s} {'M/U':>6s}")
    for (field, name), s in stacks.items():
        snr_u, snr_m = s['U']['det']['snr'], s['M']['det']['snr']
        lines.append(f"  {name + ', ' + field:28s} {snr_u:7.2f} {snr_m:7.2f} {snr_m / snr_u:6.3f}")
    text = '\n'.join(lines)
    print(text)
    with open(os.path.join(fig_dir, 'mstarweight_snr.txt'), 'w') as f:
        f.write(text + '\n')

    # LRG z2: NGC, SGC and full sky, unweighted and weighted
    fig, axes = plt.subplots(1, 2, figsize=(11, 4.2), sharey=True)
    for ax, k in zip(axes, ('U', 'M')):
        for j, (field, c) in enumerate((('full', 'k'), ('NGC', 'C0'), ('SGC', 'C3'))):
            s = stacks[(field, 'DESIY3_LRG_z2')][k]
            ax.errorbar(s['r'] + 0.05 * (j - 1), s['d'], s['sig'], fmt='o-', color=c, ms=4, capsize=2,
                        label=f"{field} (SNR {s['det']['snr']:.1f})")
        ax.axhline(0, color='gray', lw=0.8)
        ax.set_title(f"LRG z2, {'unweighted' if k == 'U' else r'weighted by $M_*^{2/3}$'}")
        ax.set_xlabel(r'$R$ [arcmin]')
        ax.legend(fontsize=9)
    axes[0].set_ylabel(r'kSZ stack, v2 filter (DSigma_lensing) [map units]')
    fig.tight_layout()
    fig.savefig(os.path.join(fig_dir, 'mstarweight_lrg_z2_caps.png'), dpi=150)

    # every full-sky sample: unweighted and weighted
    fig, axes = plt.subplots(1, len(FULL), figsize=(3.3 * len(FULL), 3.6))
    for ax, name in zip(axes, FULL):
        for j, (k, c) in enumerate((('U', 'k'), ('M', 'C1'))):
            s = stacks[('full', name)][k]
            ax.errorbar(s['r'] + 0.05 * (2 * j - 1), s['d'], s['sig'], fmt='o-', color=c, ms=3, capsize=2,
                        label=f"{'U' if k == 'U' else 'M*^(2/3)'}: SNR {s['det']['snr']:.1f}")
        ax.axhline(0, color='gray', lw=0.8)
        ax.set_title(name.replace('DESIY3_', ''))
        ax.set_xlabel(r'$R$ [arcmin]')
        ax.legend(fontsize=8)
    axes[0].set_ylabel('kSZ stack, v2 filter [map units]')
    fig.tight_layout()
    fig.savefig(os.path.join(fig_dir, 'mstarweight_full.png'), dpi=150)
    print(f"\nwrote {fig_dir}/mstarweight_snr.txt, mstarweight_lrg_z2_caps.png, mstarweight_full.png")


if __name__ == '__main__':
    sys.exit(main())
