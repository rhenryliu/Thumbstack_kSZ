"""Diagnostics of the v2-filter random-position nulls (one-off check, kept as the record behind docs/completeness_weighted_ksz_results.md), after the LRG z2
unweighted random stack failed (PTE 0.007, every aperture negative):
  1. the mean-velocity leakage of the unweighted random stack, U - U re-centred after the mask, and
     the mean filtered temperature at the random positions;
  2. a delete-one jackknife over N_PATCH contiguous sky patches (k-means, ksz_reconstruct.assign_patches),
     against ThumbStack's object bootstrap: randoms oversample the map, so neighbouring apertures
     overlap and the object bootstrap may underestimate the error;
  3. chi2 and PTE of U, U re-centred and W with the jackknife covariance (Hartlap).

Usage: python random_diagnostic.py LRG_z1 LRG_z2 ... [BGS]   (compute node; needs sklearn, e.g. the
cosmodesi_dr1 venv)
"""
import sys

import numpy as np
import pandas as pd
from scipy import stats

sys.path.insert(0, '/global/u2/r/rhliu/projects/DESIxHSC-Lensing/src')
from ksz_reconstruct import assign_patches, catalogue_mask, load_per_object_products, velocities_from_catalogue  # noqa: E402

W = '/pscratch/sd/r/rhliu/projects/Weak_lensing/'
R = W + 'ksz_measurements/ACTxDESI/spec_Y3_DSigma_lensing_compweight_random/full_no_src_with_cluster_mask_zeroV/'
CSV = {f'LRG_z{b}': (W + f'desi/spec_Y3_random/LRG_catalogues/full_zbin{b}_no_src_with_cluster_mask.csv', f'DESIY3_LRG_z{b}')
       for b in (1, 2, 3, 4)}
CSV['BGS'] = (W + 'desi/spec_Y3_random/BGS_catalogues/full_catalog_Y3_no_src_with_cluster_mask.csv', 'DESIY3_BGS')
FT = 'DSigma_lensing'
N_PATCH = 200


def ustack(t, v):
    return np.std(v) * np.sum(v[:, None] * t, axis=0) / np.sum(v ** 2)


def wstack(t, v, w):
    vp = v - np.sum(w * v) / np.sum(w)
    return np.sqrt(np.sum(w * vp ** 2) / np.sum(w)) * np.sum((w * vp)[:, None] * t, axis=0) / np.sum(w * vp ** 2)


def chi2(d, cov, n):
    c = float(d @ np.linalg.inv(cov) @ d) * (n - len(d) - 2) / (n - 1)
    return c, stats.chi2.sf(c, len(d))


def row(label, x, f='{:+.3f}'):
    print(f"  {label:44s}" + ' '.join(f.format(a) for a in x))


for name in sys.argv[1:]:
    path, sub = CSV[name]
    cat = pd.read_csv(path, usecols=['RA', 'DEC', 'VEL_LOS_RENORM', 'WEIGHT_COMP', 'WEIGHT_ZFAIL'])
    fm, fk = load_per_object_products(R + sub + '/', FT)
    assert len(fm) == len(cat), (len(fm), len(cat))
    mask = catalogue_mask(fm, fk)
    v_all = velocities_from_catalogue(cat['VEL_LOS_RENORM'].values)
    t, v = fm[mask], v_all[mask]
    w = (cat['WEIGHT_COMP'].values * cat['WEIGHT_ZFAIL'].values)[mask]
    print(f"\n=== {name} randoms, {FT}: {mask.sum()} of {len(mask)} stacked")
    U, Urc, Wt = ustack(t, v), ustack(t, v - v.mean()), wstack(t, v, w)
    ts = np.loadtxt(R + sub + f'/{FT}_ksz_uniformweight_measured.txt')
    sb = np.sqrt(np.diag(np.loadtxt(R + sub + f'/cov_{FT}_ksz_uniformweight_bootstrap.txt')))
    print(f"  offline U vs ThumbStack: max rel diff {np.max(np.abs(U / ts[:, 1] - 1)):.1e}; "
          f"mean velocity {v.mean() / v.std():+.4f} sigma_v")
    row('mean filtered temperature [map units]', t.mean(0), '{:+.5f}')
    row('U / bootstrap error', U / sb)
    row('leakage U - U re-centred / bootstrap error', (U - Urc) / sb)
    row('U re-centred / bootstrap error', Urc / sb)
    row('W / bootstrap error (unweighted)', Wt / sb)
    p = assign_patches(cat['RA'].values[mask], cat['DEC'].values[mask], n_patches=N_PATCH)
    jk = {k: [] for k in ('U', 'Urc', 'W')}
    for i in range(N_PATCH):
        keep = p != i
        tt, vv, ww = t[keep], v[keep], w[keep]
        jk['U'].append(ustack(tt, vv))
        jk['Urc'].append(ustack(tt, vv - vv.mean()))
        jk['W'].append(wstack(tt, vv, ww))
    for k, x in (('U', U), ('Urc', Urc), ('W', Wt)):
        a = np.array(jk[k])
        cov = (N_PATCH - 1) / N_PATCH * np.einsum('ij,ik->jk', a - a.mean(0), a - a.mean(0))
        c, pte = chi2(x, cov, N_PATCH)
        row(f'{k}: jackknife / unweighted bootstrap error', np.sqrt(np.diag(cov)) / sb, '{:.2f}')
        print(f"  {k}: chi2 with the jackknife covariance {c:.1f} / {len(x)}, PTE {pte:.3f}")
