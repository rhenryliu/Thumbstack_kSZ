"""Known answer for ThumbStack's 'ksz_mstarweight' (one-off check, kept as the record of the configs/mstarweight/ run).

Recomputes the stacks offline from the per-object filtering the run read (filtering_dir), with the catalogue
mask of ThumbStack.catalogMask (point-source cut; objects with LOGMSTAR <= 0 dropped; outlier cut on the kept
objects), and compares:
  - ThumbStack ksz_mstarweight   vs offline weighted, w = M*^(2/3)  (goal: <= 1e-12 relative)
  - ThumbStack ksz_uniformweight vs offline unweighted              (goal: <= 1e-12 relative)
Then splits the difference between this run's unweighted stack and the earlier one of all objects into
the objects dropped and the outlier threshold, with and without the mean-velocity leakage (D4: the
unweighted estimator does not re-centre the velocities on the stacked objects).

Usage: python mstarweight_known_answer.py [LRG_z1 ... BGS LRG_z2_NGC LRG_z2_SGC]  (compute node;
cosmodesi_dr1 venv, Python >= 3.10 for DESIxHSC-Lensing src/ksz_reconstruct.py)
Recorded run: all samples, job 61616995 (log /pscratch/sd/r/rhliu/logs/mstarweight-ka-20261119-161016.log).
In it the outlier cut removes the same objects with and without the LOGMSTAR objects ("outlier threshold
alone: 0 out, 0 in" for every sample), so the run checks the estimators but cannot tell the order of the
cuts apart, and the "objects dropped" line equals the total difference by construction. The "mean T of
dropped / mean T of all" ratio is ill-conditioned where the mean temperature is near zero.
"""
import argparse
import sys

import numpy as np
import pandas as pd
from scipy import optimize, special

sys.path.insert(0, '/global/u2/r/rhliu/projects/DESIxHSC-Lensing/src')
from ksz_reconstruct import ksz_delta_sigma, load_per_object_products, velocities_from_catalogue  # noqa: E402

W = '/pscratch/sd/r/rhliu/projects/Weak_lensing/'
K = W + 'ksz_measurements/ACTxDESI/'
FT = 'DSigma_lensing'
CUT = 'no_src_with_cluster_mask'
NEW = K + 'spec_Y3_DSigma_lensing_mstarweight/'
POWER = 2. / 3.
SAMPLES = {f'LRG_z{b}': dict(cat=W + f'desi/spec_Y3/LRG_catalogues_zeroV/full_zbin{b}_{CUT}.csv',
                             filt=K + f'spec_Y3_DSigma_lensing/full_{CUT}_zeroV/DESIY3_LRG_z{b}/',
                             sub=f'full_{CUT}_zeroV/DESIY3_LRG_z{b}/') for b in (1, 2, 3, 4)}
SAMPLES['BGS'] = dict(cat=W + f'desi/spec_Y3/BGS_CIGALE/BGS_catalogues/full_catalog_Y3_{CUT}.csv',
                      filt=K + f'spec_Y3_DSigma_lensing/CIGALE/full_{CUT}_zeroV/DESIY3_BGS/',
                      sub=f'CIGALE/full_{CUT}_zeroV/DESIY3_BGS/')
for cap in ('NGC', 'SGC'):
    SAMPLES[f'LRG_z2_{cap}'] = dict(cat=W + f'desi/spec_Y3/LRG_catalogues_zeroV/{cap}_zbin2_{CUT}.csv',
                                    filt=K + f'spec_Y3_DSigma_lensing_compweight/{cap}_{CUT}_zeroV/DESIY3_LRG_z2/',
                                    sub=f'{cap}_{CUT}_zeroV/DESIY3_LRG_z2/')


def thumbstack_mask(t, fmask, keep):
    """ThumbStack.catalogMask: point-source cut and *keep*, then the outlier cut on those objects."""
    mask = (np.abs(fmask[:, -1]) < 1.) & keep
    n_obj = int(mask.sum())
    target = special.erfc(5. / np.sqrt(2.))
    n_cut = optimize.brentq(lambda n: n_obj * special.erfc(n / np.sqrt(2.)) - target, 0., 1.e2)
    sigmas = np.std(t[mask, :], axis=0)
    return mask & np.all(np.abs(t) <= n_cut * sigmas[None, :], axis=1)


def wstack(t, v, w):
    """Weighted estimator from its definition: sigma_w * sum w v' T / sum w v'^2, v' = v - <v>_w."""
    vp = v - np.sum(w * v) / np.sum(w)
    sw = np.sqrt(np.sum(w * vp**2) / np.sum(w))
    return sw * np.sum((w * vp)[:, None] * t, axis=0) / np.sum(w * vp**2)


def fmt(x, f='{:+.3f}'):
    return ' '.join(f.format(v) for v in x)


def main():
    p = argparse.ArgumentParser()
    p.add_argument('samples', nargs='*', default=list(SAMPLES))
    for name in p.parse_args().samples:
        s = SAMPLES[name]
        cat = pd.read_csv(s['cat'], usecols=['VEL_LOS_RENORM', 'LOGMSTAR'])
        t, fmask = load_per_object_products(s['filt'], FT)
        if len(t) != len(cat):
            raise ValueError(f"{name}: {len(t)} filtered objects, {len(cat)} catalogue rows")
        v = velocities_from_catalogue(cat['VEL_LOS_RENORM'].values)
        logm = cat['LOGMSTAR'].values.astype(float)
        good = np.isfinite(logm) & (logm > 0.)
        new = thumbstack_mask(t, fmask, good)
        old = thumbstack_mask(t, fmask, np.ones(len(t), bool))
        w = 10.**(POWER * logm[new])
        ref_m, ref_u = wstack(t[new], v[new], w), ksz_delta_sigma(t[new], v[new])
        d = NEW + s['sub']
        ts_m = np.loadtxt(d + f'{FT}_ksz_mstarweight_measured.txt')[:, 1]
        ts_u = np.loadtxt(d + f'{FT}_ksz_uniformweight_measured.txt')[:, 1]
        sig = np.sqrt(np.diag(np.loadtxt(d + f'cov_{FT}_ksz_uniformweight_bootstrap.txt')))
        earlier = np.loadtxt(s['filt'] + f'{FT}_ksz_uniformweight_measured.txt')[:, 1]
        print(f"\n=== {name}: {new.sum()} stacked here, {old.sum()} before; dropped for LOGMSTAR <= 0: "
              f"{(old & ~good).sum()}; outlier threshold alone: {(old & good & ~new).sum()} out, "
              f"{(new & ~old).sum()} in")
        print(f"  ksz_mstarweight vs offline weighted: max rel diff {np.max(np.abs(ts_m / ref_m - 1)):.2e}")
        print(f"  ksz_uniformweight vs offline unweighted: max rel diff {np.max(np.abs(ts_u / ref_u - 1)):.2e}")
        print(f"  earlier U vs offline unweighted of all objects: max rel diff "
              f"{np.max(np.abs(earlier / ksz_delta_sigma(t[old], v[old]) - 1)):.2e}")
        # the difference to the earlier stack, in units of this run's unweighted bootstrap sigma
        same_thr = old & good  # the earlier mask with the objects without mass removed
        print("  (U - earlier U) / sigma                  " + fmt((ts_u - earlier) / sig))
        print("    objects dropped, earlier threshold     " + fmt((ksz_delta_sigma(t[same_thr], v[same_thr]) - earlier) / sig))
        print("    both re-centred (no D4 leakage)        " + fmt((ksz_delta_sigma(t[new], v[new], rezero=True)
                                                                   - ksz_delta_sigma(t[old], v[old], rezero=True)) / sig))
        sv = v[old].std()
        print(f"  mean velocity of the stacked objects: before {v[old].mean() / sv:+.4f} sigma_v, here "
              f"{v[new].mean() / sv:+.4f} sigma_v; of the dropped objects {v[old & ~good].mean() / sv:+.3f} sigma_v")
        drop = old & ~good
        print("  dropped objects: mean T / mean T of all  " + fmt(t[drop].mean(axis=0) / t[old].mean(axis=0), '{:.2f}'))


if __name__ == '__main__':
    main()
