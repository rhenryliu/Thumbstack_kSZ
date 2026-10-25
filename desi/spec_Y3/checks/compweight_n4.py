"""N4, weight controls for the completeness-weighted kSZ stacks (one-off check, kept as the record behind docs/completeness_weighted_ksz_results.md).

From the production per-object filtering (full sky), for each sample and filter:
  - weights shuffled among the stacked galaxies, velocities fixed (N_PERM permutations): the GLS
    amplitude W_perm / U must scatter around 1; the observed W / U is placed in that distribution;
  - leakage without re-centring the velocities on the weighted mean, (W0 - W) / sigma_U;
  - the weighted and unweighted mean temperature per aperture;
  - D4: the mean velocity of the stacked objects, and the offset it predicts for the unweighted
    within-sample shuffle null, (<v>/sigma_v) <T> / sigma_U per aperture, to compare with N1 (b);
  - an independent within-sample shuffle null of the weighted stack (N_SHUF2 realisations; independent
    of ThumbStack's shuffles, but drawn from the same rng as the weight permutations, so the results
    depend on the argument order), to tell a fluctuation of N1 (b) from a bias.
sigma_U and the GLS covariance are ThumbStack's unweighted bootstrap covariance (compweight folders).

Usage: python compweight_n4.py [LRG_z1 ... BGS] [--filters DSigma_lensing DSigma]  (compute node)
Recorded run (job 60178685): compweight_n4.py LRG_z1 LRG_z2 LRG_z3 LRG_z4 BGS --filters DSigma_lensing
"""
import argparse

import numpy as np

from compweight_known_answer import NEW, PROD, SAMPLES, catalogue_mask, load_per_object_products, wstack
from compweight_known_answer import fitsio, pd, SIDE, velocities_from_catalogue

N_PERM = 200
N_SHUF2 = 2000
rng = np.random.default_rng(4)


def ustack(t, v):
    """ThumbStack's ksz_uniformweight (no re-centring after the mask)."""
    return np.std(v) * np.sum(v[:, None] * t, axis=0) / np.sum(v ** 2)


def w0stack(t, v, w):
    """Weighted, velocities not re-centred on the weighted mean."""
    vm = np.sum(w * v) / np.sum(w)
    sw = np.sqrt(np.sum(w * (v - vm) ** 2) / np.sum(w))
    return sw * np.sum((w * v)[:, None] * t, axis=0) / np.sum(w * v ** 2)


def row(label, x, f='{:+.3f}'):
    print(f"  {label:46s}" + ' '.join(f.format(a) for a in x))


def main():
    p = argparse.ArgumentParser()
    p.add_argument('samples', nargs='*', default=list(SAMPLES))
    p.add_argument('--filters', nargs='+', default=['DSigma_lensing'])  # the v1 filter only when asked
    a = p.parse_args()
    for name in a.samples:
        s = SAMPLES[name]
        if s.get('side'):
            cat = pd.read_csv(s['cat'], usecols=['VEL_LOS_RENORM'])
            side = fitsio.read(SIDE, columns=['WEIGHT_COMP', 'WEIGHT_ZFAIL'])
            w_all = side['WEIGHT_COMP'].astype(float) * side['WEIGHT_ZFAIL'].astype(float)
        else:
            cat = pd.read_csv(s['cat'], usecols=['VEL_LOS_RENORM', 'WEIGHT_COMP', 'WEIGHT_ZFAIL'])
            w_all = cat['WEIGHT_COMP'].values.astype(float) * cat['WEIGHT_ZFAIL'].values.astype(float)
        v_all = velocities_from_catalogue(cat['VEL_LOS_RENORM'].values)
        for ft in a.filters:
            print(f"\n=== {name}, {ft}")
            fm, fk = load_per_object_products(PROD[ft] + s['sub'], ft)
            mask = catalogue_mask(fm, fk)
            t, v, w = fm[mask], v_all[mask], w_all[mask]
            try:
                cu = np.loadtxt(NEW[ft] + s['sub'] + f'cov_{ft}_ksz_uniformweight_bootstrap.txt')
            except OSError as e:
                print(f"  no ThumbStack covariance yet: {e}")
                continue
            su, ci = np.sqrt(np.diag(cu)), np.linalg.inv(cu)
            U, Wr = ustack(t, v), wstack(t, v, w)
            amp = lambda x: (U @ ci @ x) / (U @ ci @ U)  # noqa: E731
            perm = np.array([amp(wstack(t, v, rng.permutation(w))) for _ in range(N_PERM)])
            print(f"  amplitude W/U {amp(Wr):.3f}; weights shuffled among galaxies: {perm.mean():.3f} +- {perm.std():.3f} "
                  f"({N_PERM}); observed is {(amp(Wr) - perm.mean()) / perm.std():+.1f} sigma from it")
            row('leakage without re-centring (W0 - W) / sigma_U', (w0stack(t, v, w) - Wr) / su)
            row('mean temperature, unweighted [map units]', t.mean(0), '{:+.4f}')
            row('mean temperature, weighted', np.sum(w[:, None] * t, 0) / w.sum(), '{:+.4f}')
            vb = v.mean() / v.std()
            print(f"  D4: mean velocity of the stacked objects {vb:+.4f} sigma_v (weighted {np.sum(w * v) / w.sum() / v.std():+.4f})")
            row('D4: predicted unweighted within-sample null / sigma_U', vb * t.mean(0) / su)
            row('D4: same, ThumbStack U offset vs U re-centred / sigma_U', (U - ustack(t, v - v.mean())) / su)
            sh = np.array([wstack(t, rng.permutation(v), w) for _ in range(N_SHUF2)])
            m, c = sh.mean(0), np.cov(sh, rowvar=False) / N_SHUF2
            chi2 = m @ np.linalg.inv(c) @ m * (N_SHUF2 - len(m) - 2) / (N_SHUF2 - 1)
            from scipy import stats
            print(f"  independent within-sample shuffle, weighted ({N_SHUF2}): chi2 {chi2:.1f} / {len(m)}, "
                  f"PTE {stats.chi2.sf(chi2, len(m)):.3f}")
            row('  its mean / sigma_W', m / np.sqrt(np.diag(np.loadtxt(NEW[ft] + s['sub'] + f'cov_{ft}_ksz_compweight_bootstrap.txt'))))


if __name__ == '__main__':
    main()
