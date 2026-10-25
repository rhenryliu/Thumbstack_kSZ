"""Null tests of the completeness-weighted kSZ stacks (estimator 'ksz_compweight'), side by side with the
unweighted ones ('ksz_uniformweight'), from the outputs of configs/compweight/.

Pass criteria, fixed before the runs: PTE >= 0.05 for each sample and filter, and, for the means over
shuffled-velocity realisations, no aperture above 0.3 sigma, sigma being the bootstrap error of the
measured stack (same sample, filter and estimator).

  Stacks    weighted / unweighted per aperture (paired bootstrap error, same resamples) and in
            amplitude (GLS with the unweighted bootstrap covariance), and against PIP weights.
  N1        shuffled velocities: (a) paper style, whole catalogue shuffled before the mask
            (make_stacks_shuffle.py, n_seeds realisations); (b) within the stacked objects (ThumbStack's
            vshuffle). chi2 of the mean against zero with covariance C/N (Hartlap, N realisations), and
            the per-realisation chi2 against chi2(nRAp).
  N2        random positions: chi2 of the stack against zero with its bootstrap covariance (Hartlap).
  N3        NGC against SGC: chi2 of the difference with C_NGC + C_SGC, each Hartlap corrected
            (as DESIxHSC-Lensing src/snr.py consistency_test, paper Figure 25).

Usage: python compweight_nulls.py [--filterType DSigma_lensing]
"""
import argparse
import os

import numpy as np
from scipy import stats

K = '/pscratch/sd/r/rhliu/projects/Weak_lensing/ksz_measurements/ACTxDESI/'
ROOT = {'DSigma_lensing': K + 'spec_Y3_DSigma_lensing_compweight', 'DSigma': K + 'spec_Y3_DSigma_compweight'}
CUT = 'no_src_with_cluster_mask_zeroV'
LRG = [f'DESIY3_LRG_z{b}' for b in (1, 2, 3, 4)]
EST = ('ksz_uniformweight', 'ksz_compweight')
N_BOOT = 10000
PTE_MIN, MEAN_MAX = 0.05, 0.3


def sample_dir(root, name, field='full', cigale=True):
    """Output folder of one stack; the full-sky BGS stacks use the CIGALE catalogue."""
    sub = 'CIGALE/' if (name == 'DESIY3_BGS' and cigale) else ''
    return f'{root}/{sub}{field}_{CUT}/{name}/'


def load(path):
    return np.loadtxt(path) if os.path.exists(path) else None


def hartlap(n, p):
    return (n - p - 2) / (n - 1)


def chi2_pte(d, cov, n_resample):
    """chi2 of d against zero with the Hartlap-corrected inverse of cov, and its PTE."""
    chi2 = float(d @ np.linalg.inv(cov) @ d) * hartlap(n_resample, len(d))
    return chi2, float(stats.chi2.sf(chi2, df=len(d)))


def flag(ok):
    return 'pass' if ok else 'FAIL'


def ratio(d, ft):
    """Weighted / unweighted stack of one folder: per aperture, with the linearised paired-bootstrap
    error, and as a GLS amplitude with the unweighted covariance; None if a file is missing."""
    prof = {e: load(d + f'{ft}_{e}_measured.txt') for e in EST}
    cov = {e: load(d + f'cov_{ft}_{e}_bootstrap.txt') for e in EST}
    boot = {e: load(d + f'{ft}_{e}_bootstrap_samples.txt') for e in EST}
    if any(x is None for x in list(prof.values()) + list(cov.values()) + list(boot.values())):
        return None
    u, w = prof['ksz_uniformweight'][:, 1], prof['ksz_compweight'][:, 1]
    bu, bw = boot['ksz_uniformweight'], boot['ksz_compweight']
    ci = np.linalg.inv(cov['ksz_uniformweight'])
    amps = np.einsum('ij,jk,ik->i', bu, ci, bw) / np.einsum('ij,jk,ik->i', bu, ci, bu)
    return dict(r=w / u, err=np.std((bw - (w / u) * bu) / u, axis=0), amp=(u @ ci @ w) / (u @ ci @ u),
                amp_err=np.std(amps), sig={e: np.sqrt(np.diag(cov[e])) for e in EST})


def stacks(ft, root):
    print(f"\n## Stacks, {ft}: weighted / unweighted")
    sig = {}
    for name in LRG + ['DESIY3_BGS']:
        d = sample_dir(root, name)
        prof = {e: load(d + f'{ft}_{e}_measured.txt') for e in EST}
        cov = {e: load(d + f'cov_{ft}_{e}_bootstrap.txt') for e in EST}
        boot = {e: load(d + f'{ft}_{e}_bootstrap_samples.txt') for e in EST}
        if any(prof[e] is None or cov[e] is None for e in EST):
            print(f"  {name}: missing")
            continue
        sig[name] = {e: np.sqrt(np.diag(cov[e])) for e in EST}
        u, w = prof['ksz_uniformweight'][:, 1], prof['ksz_compweight'][:, 1]
        ci = np.linalg.inv(cov['ksz_uniformweight'])
        amp = (u @ ci @ w) / (u @ ci @ u)
        line = f"  {name}: W/U " + ' '.join(f'{x:.3f}' for x in w / u)
        if boot['ksz_uniformweight'] is not None and boot['ksz_compweight'] is not None:
            bu, bw = boot['ksz_uniformweight'], boot['ksz_compweight']
            amps = np.einsum('ij,jk,ik->i', bu, ci, bw) / np.einsum('ij,jk,ik->i', bu, ci, bu)
            # linearised (delta-method) paired error: the ratio of resamples is heavy tailed where U is noisy
            err = np.std((bw - (w / u) * bu) / u, axis=0)
            line += "\n      paired-bootstrap error " + ' '.join(f'{x:.3f}' for x in err)
            line += f"\n      amplitude W/U (GLS, unweighted covariance) {amp:.3f} +- {np.std(amps):.3f}"
        else:
            line += f"\n      amplitude W/U (GLS, unweighted covariance) {amp:.3f}"
        line += "\n      bootstrap error W/U " + ' '.join(f'{x:.3f}' for x in sig[name]['ksz_compweight'] / sig[name]['ksz_uniformweight'])
        print(line)
        # the PIP run drops the objects without a PIP weight, from both estimators
        pip = ratio(sample_dir(root + '_pip', name), ft)
        if pip is not None:
            print("      PIP W / PIP U " + ' '.join(f'{x:.3f}' for x in pip['r'])
                  + f"\n      PIP amplitude W/U {pip['amp']:.3f} +- {pip['amp_err']:.3f}")
    return sig


def n1(ft, root, sig):
    print(f"\n## N1, shuffled velocities, {ft}")
    for name in LRG + ['DESIY3_BGS']:
        if name not in sig:
            continue
        for e in EST:
            # (a) paper style: the shuffle driver writes the mean over seeds and C/N
            d = sample_dir(root + '_shuffle', name)
            mean, covm, samples = load(d + f'{ft}_{e}_measured.txt'), load(d + f'cov_{ft}_{e}_bootstrap.txt'), \
                load(d + f'{ft}_{e}_shuffle_samples.txt')
            if mean is not None and covm is not None and samples is not None:
                n = len(samples)
                report(f"{name} {e} (a) whole-catalogue shuffle, {n} seeds", mean[:, 1], covm, n, samples, sig[name][e])
            else:
                print(f"  {name} {e} (a): missing")
            # (b) within the stacked objects
            d = sample_dir(root, name)
            mean, cov1, samples = load(d + f'{ft}_{e}_vshufflemean.txt'), load(d + f'cov_{ft}_{e}_vshuffle.txt'), \
                load(d + f'{ft}_{e}_vshuffle_samples.txt')
            if mean is not None and cov1 is not None and samples is not None:
                n = len(samples)
                report(f"{name} {e} (b) within-sample shuffle, {n}", mean[:, 1], cov1 / n, n, samples, sig[name][e])
            else:
                print(f"  {name} {e} (b): missing")


def report(label, mean, cov_mean, n, samples, sigma):
    chi2, pte = chi2_pte(mean, cov_mean, n)
    worst = np.max(np.abs(mean) / sigma)
    s = f"  {label}: [{flag(pte >= PTE_MIN and worst <= MEAN_MAX)}] chi2 {chi2:.1f} / {len(mean)}, PTE {pte:.3f}; " \
        f"max |mean| {worst:.3f} sigma"
    if samples is not None:
        # each realisation against zero, with the scatter of the realisations as covariance (no Hartlap
        # factor: each realisation is in the covariance it is tested against, which already biases low)
        c1 = np.cov(samples, rowvar=False)
        per = np.einsum('ij,jk,ik->i', samples, np.linalg.inv(c1), samples)
        s += f"; per-realisation chi2: mean {per.mean():.2f} (expect {samples.shape[1]}), " \
             f"KS p {stats.kstest(per, 'chi2', args=(samples.shape[1],)).pvalue:.3f}"
    print(s + "\n      mean / sigma " + ' '.join(f'{x:+.3f}' for x in mean / sigma))


def n2(ft, root):
    print(f"\n## N2, random positions, {ft} (DSigma_lensing randoms only)")
    for wlab, r in (('WEIGHT_COMP x WEIGHT_ZFAIL', root + '_random'), ('WEIGHT_COMP', root + '_random_comp')):
        for name in LRG + ['DESIY3_BGS']:
            d = sample_dir(r, name, cigale=False)
            for e in EST:
                if e == 'ksz_uniformweight' and wlab == 'WEIGHT_COMP':
                    continue  # same as above
                prof, cov = load(d + f'{ft}_{e}_measured.txt'), load(d + f'cov_{ft}_{e}_bootstrap.txt')
                if prof is None or cov is None:
                    continue
                chi2, pte = chi2_pte(prof[:, 1], cov, N_BOOT)
                print(f"  {name} {e} [{wlab}]: chi2 {chi2:.1f} / {len(prof)}, PTE {pte:.3f} [{flag(pte >= PTE_MIN)}]"
                      "\n      stack / own error " + ' '.join(f'{x:+.3f}' for x in prof[:, 1] / np.sqrt(np.diag(cov))))


def n3(ft, root):
    print(f"\n## N3, NGC against SGC, {ft} (BGS: non-CIGALE catalogue, as Figure 25)")
    for name in LRG + ['DESIY3_BGS']:
        for e in EST:
            a, b = (sample_dir(root, name, field=f, cigale=False) for f in ('NGC', 'SGC'))
            pa, ca, pb, cb = (load(a + f'{ft}_{e}_measured.txt'), load(a + f'cov_{ft}_{e}_bootstrap.txt'),
                              load(b + f'{ft}_{e}_measured.txt'), load(b + f'cov_{ft}_{e}_bootstrap.txt'))
            if any(x is None for x in (pa, ca, pb, cb)):
                continue
            h = hartlap(N_BOOT, len(pa))
            delta = pa[:, 1] - pb[:, 1]
            chi2 = float(delta @ np.linalg.inv(ca / h + cb / h) @ delta)
            pte = float(stats.chi2.sf(chi2, df=len(delta)))
            print(f"  {name} {e}: chi2 {chi2:.1f} / {len(delta)}, PTE {pte:.3f} [{flag(pte >= PTE_MIN)}]")


def main():
    p = argparse.ArgumentParser(description=__doc__.split('\n')[0])
    p.add_argument('--filterType', default='DSigma_lensing', choices=list(ROOT))
    ft = p.parse_args().filterType
    root = ROOT[ft]
    sig = stacks(ft, root)
    n1(ft, root, sig)
    n2(ft, root)
    n3(ft, root)


if __name__ == '__main__':
    main()
