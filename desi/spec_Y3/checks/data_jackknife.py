"""Bootstrap against jackknife errors for the v2-filter (DSigma_lensing) DATA stacks (one-off check, 2026-10-25;
output in ksz_measurements/ACTxDESI/restack_checks/data_jackknife_20261025.txt, job 60253508, run
from a scratch copy identical to this file but for these two docstring lines). Each jackknife error
ratio carries about sqrt(1/(n_patch - 1)) noise (about 5% for 200 patches); the jackknife S/N
includes the Hartlap factor (0.975 for 200 patches).

The random-position nulls showed ThumbStack's object bootstrap underestimating the errors by up to
x1.35 at large apertures (desi/spec_Y3/checks/random_diagnostic.py). This repeats that comparison
on the data stacks, unweighted (ksz_uniformweight, production) and weighted (ksz_compweight):
  1. offline U and W from the production per-object filtering, against ThumbStack's measured files
     (known answer);
  2. delete-one jackknife over N_PATCH contiguous sky patches (k-means, ksz_reconstruct.assign_patches),
     for N_PATCH = 100 and 200, against the bootstrap covariances ThumbStack saved;
  3. S/N (sqrt chi2 against zero, Hartlap) with each covariance.
Writes nothing; print only.

Usage: python data_jackknife.py [LRG_z1 ... BGS]   (compute node)
"""
import sys

import fitsio
import numpy as np
import pandas as pd

sys.path.insert(0, '/global/u2/r/rhliu/projects/DESIxHSC-Lensing/src')
from ksz_reconstruct import assign_patches, catalogue_mask, load_per_object_products, velocities_from_catalogue  # noqa: E402

W = '/pscratch/sd/r/rhliu/projects/Weak_lensing/'
K = W + 'ksz_measurements/ACTxDESI/'
PROD = K + 'spec_Y3_DSigma_lensing/'
NEW = K + 'spec_Y3_DSigma_lensing_compweight/'
SIDE = W + 'desi/spec_Y3/lensing/CIGALE/lens_weights/desi_BGS_bin0.fits'
SAMPLES = {f'LRG_z{b}': dict(cat=W + f'desi/spec_Y3/LRG_catalogues_zeroV/full_zbin{b}_no_src_with_cluster_mask.csv',
                             sub=f'full_no_src_with_cluster_mask_zeroV/DESIY3_LRG_z{b}/') for b in (1, 2, 3, 4)}
SAMPLES['BGS'] = dict(cat=W + 'desi/spec_Y3/BGS_CIGALE/BGS_catalogues/full_catalog_Y3_no_src_with_cluster_mask.csv',
                      sub='CIGALE/full_no_src_with_cluster_mask_zeroV/DESIY3_BGS/', side=True)
FT = 'DSigma_lensing'


def ustack(t, v):
    """ThumbStack's ksz_uniformweight (no re-centring after the mask)."""
    return np.std(v) * np.sum(v[:, None] * t, axis=0) / np.sum(v ** 2)


def wstack(t, v, w):
    """ThumbStack's ksz_compweight (re-centred on the weighted mean)."""
    vp = v - np.sum(w * v) / np.sum(w)
    return np.sqrt(np.sum(w * vp ** 2) / np.sum(w)) * np.sum((w * vp)[:, None] * t, axis=0) / np.sum(w * vp ** 2)


def snr(x, cov, n):
    """sqrt(chi2) against zero with the Hartlap factor for n realisations."""
    return np.sqrt(float(x @ np.linalg.inv(cov) @ x) * (n - len(x) - 2) / (n - 1))


def row(label, x, f='{:.2f}'):
    print(f"  {label:46s}" + ' '.join(f.format(a) for a in x))


for name in (sys.argv[1:] or list(SAMPLES)):
    s = SAMPLES[name]
    if s.get('side'):
        cat = pd.read_csv(s['cat'], usecols=['RA', 'DEC', 'Z', 'VEL_LOS_RENORM'])
        side = fitsio.read(SIDE)
        for c in ('RA', 'DEC', 'Z'):
            assert np.abs(side[c].astype(float) - cat[c].values).max() < 1e-9, c
        w_all = side['WEIGHT_COMP'].astype(float) * side['WEIGHT_ZFAIL'].astype(float)
    else:
        cat = pd.read_csv(s['cat'], usecols=['RA', 'DEC', 'VEL_LOS_RENORM', 'WEIGHT_COMP', 'WEIGHT_ZFAIL'])
        w_all = cat['WEIGHT_COMP'].values.astype(float) * cat['WEIGHT_ZFAIL'].values.astype(float)
    fm, fk = load_per_object_products(PROD + s['sub'], FT)
    assert len(fm) == len(cat), (len(fm), len(cat))
    mask = catalogue_mask(fm, fk)
    v_all = velocities_from_catalogue(cat['VEL_LOS_RENORM'].values)
    t, v, w = fm[mask], v_all[mask], w_all[mask]
    print(f"\n=== {name}, {FT}: {mask.sum()} of {len(mask)} stacked")
    U, Wt = ustack(t, v), wstack(t, v, w)
    ts_u = np.loadtxt(PROD + s['sub'] + f'{FT}_ksz_uniformweight_measured.txt')[:, 1]
    ts_w = np.loadtxt(NEW + s['sub'] + f'{FT}_ksz_compweight_measured.txt')[:, 1]
    print(f"  known answer: offline U vs ThumbStack max rel diff {np.max(np.abs(U / ts_u - 1)):.1e}; "
          f"offline W vs ThumbStack {np.max(np.abs(Wt / ts_w - 1)):.1e}")
    cb = {'U': np.loadtxt(PROD + s['sub'] + f'cov_{FT}_ksz_uniformweight_bootstrap.txt'),
          'W': np.loadtxt(NEW + s['sub'] + f'cov_{FT}_ksz_compweight_bootstrap.txt')}
    x = {'U': U, 'W': Wt}
    for k in ('U', 'W'):
        print(f"  {k}: S/N with the bootstrap covariance {snr(x[k], cb[k], 10000):.2f}")
    for n_patch in (100, 200):
        p = assign_patches(cat['RA'].values[mask], cat['DEC'].values[mask], n_patches=n_patch)
        jk = {'U': [], 'W': []}
        for i in range(n_patch):
            keep = p != i
            jk['U'].append(ustack(t[keep], v[keep]))
            jk['W'].append(wstack(t[keep], v[keep], w[keep]))
        for k in ('U', 'W'):
            a = np.array(jk[k])
            cov = (n_patch - 1) / n_patch * np.einsum('ij,ik->jk', a - a.mean(0), a - a.mean(0))
            row(f'{k}, {n_patch} patches: jackknife / bootstrap error', np.sqrt(np.diag(cov) / np.diag(cb[k])))
            print(f"  {k}, {n_patch} patches: S/N with the jackknife covariance {snr(x[k], cov, n_patch):.2f}")
