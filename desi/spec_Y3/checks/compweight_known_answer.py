"""Known answer for ThumbStack's 'ksz_compweight' (one-off check, kept as the record behind docs/completeness_weighted_ksz_results.md).

Recomputes the weighted and unweighted stacks offline from the production per-object filtering, with
DESIxHSC-Lensing's src/ksz_reconstruct.py and the weighted estimator of
DESIxHSC-Lensing-v2.2/scratch/v2.2_ksz_weights/ksz_weight_check.py (wstack, copied below), and compares:
  - ThumbStack ksz_compweight (new folder)   vs offline weighted    (goal: <= 1e-12 relative)
  - ThumbStack ksz_uniformweight (new folder) vs production profile  (goal: byte-identical files)
Also prints the D4 diagnostic: the mean velocity of the stacked objects, unweighted and weighted.

Usage: python compweight_known_answer.py [LRG_z1 BGS ...] [--filters DSigma_lensing DSigma]  (compute node)
Recorded runs: 'LRG_z1 BGS' and 'LRG_z2 LRG_z3 LRG_z4', both with --filters DSigma_lensing.
"""
import argparse
import filecmp
import sys

import fitsio
import numpy as np
import pandas as pd

sys.path.insert(0, '/global/u2/r/rhliu/projects/DESIxHSC-Lensing/src')
from ksz_reconstruct import catalogue_mask, ksz_delta_sigma, load_per_object_products, velocities_from_catalogue  # noqa: E402

W = '/pscratch/sd/r/rhliu/projects/Weak_lensing/'
K = W + 'ksz_measurements/ACTxDESI/'
SIDE = W + 'desi/spec_Y3/lensing/CIGALE/lens_weights/desi_BGS_bin0.fits'
PROD = {'DSigma_lensing': K + 'spec_Y3_DSigma_lensing/', 'DSigma': K + 'spec_Y3/'}
NEW = {'DSigma_lensing': K + 'spec_Y3_DSigma_lensing_compweight/', 'DSigma': K + 'spec_Y3_DSigma_compweight/'}
SAMPLES = {f'LRG_z{b}': dict(cat=W + f'desi/spec_Y3/LRG_catalogues_zeroV/full_zbin{b}_no_src_with_cluster_mask.csv',
                             sub=f'full_no_src_with_cluster_mask_zeroV/DESIY3_LRG_z{b}/') for b in (1, 2, 3, 4)}
SAMPLES['BGS'] = dict(cat=W + 'desi/spec_Y3/BGS_CIGALE/BGS_catalogues/full_catalog_Y3_no_src_with_cluster_mask.csv',
                      sub='CIGALE/full_no_src_with_cluster_mask_zeroV/DESIY3_BGS/', side=True)


def wstack(t, v, w):
    """ksz_weight_check.py's weighted estimator (re-zeroed on the weighted mean)."""
    vm = np.sum(w * v) / np.sum(w)
    vp = v - vm
    sw = np.sqrt(np.sum(w * (v - vm) ** 2) / np.sum(w))
    return sw * np.sum(w[:, None] * vp[:, None] * t, axis=0) / np.sum(w * vp ** 2)


def main():
    p = argparse.ArgumentParser()
    p.add_argument('samples', nargs='*', default=list(SAMPLES))
    p.add_argument('--filters', nargs='+', default=['DSigma_lensing'])  # the v1 filter only when asked
    a = p.parse_args()
    for name in a.samples:
        s = SAMPLES[name]
        if s.get('side'):
            cat = pd.read_csv(s['cat'], usecols=['RA', 'DEC', 'Z', 'VEL_LOS_RENORM'])
            side = fitsio.read(SIDE)
            for c in ('RA', 'DEC', 'Z'):
                assert np.abs(side[c].astype(float) - cat[c].values).max() < 1e-9, c
            w_all = side['WEIGHT_COMP'].astype(float) * side['WEIGHT_ZFAIL'].astype(float)
        else:
            cat = pd.read_csv(s['cat'], usecols=['VEL_LOS_RENORM', 'WEIGHT_COMP', 'WEIGHT_ZFAIL'])
            w_all = cat['WEIGHT_COMP'].values.astype(float) * cat['WEIGHT_ZFAIL'].values.astype(float)
        v_all = velocities_from_catalogue(cat['VEL_LOS_RENORM'].values)
        for ft in a.filters:
            prod, new = PROD[ft] + s['sub'], NEW[ft] + s['sub']
            print(f"\n=== {name}, {ft}")
            fm, fk = load_per_object_products(prod, ft)
            mask = catalogue_mask(fm, fk)
            t, v, w = fm[mask], v_all[mask], w_all[mask]
            sv = v.std()
            print(f"{mask.sum()} of {len(mask)} objects stacked; mean velocity (D4): unweighted {v.mean() / sv:+.4f} sigma_v, "
                  f"weighted {np.sum(w * v) / np.sum(w) / sv:+.4f} sigma_v; mean weight {w.mean():.4f}")
            ref_w, ref_u = wstack(t, v, w), ksz_delta_sigma(t, v)
            try:
                ts_w = np.loadtxt(new + f'{ft}_ksz_compweight_measured.txt')[:, 1]
                ts_u_file = new + f'{ft}_ksz_uniformweight_measured.txt'
                ts_u = np.loadtxt(ts_u_file)[:, 1]
            except OSError as e:
                print(f"  no ThumbStack output yet: {e}")
                continue
            print(f"  ksz_compweight vs offline weighted: max rel diff {np.max(np.abs(ts_w / ref_w - 1)):.2e}")
            print(f"  ksz_uniformweight vs offline unweighted: max rel diff {np.max(np.abs(ts_u / ref_u - 1)):.2e}")
            for f in (f'{ft}_ksz_uniformweight_measured.txt', f'cov_{ft}_ksz_uniformweight_bootstrap.txt'):
                print(f"  {f}: new == production byte for byte: {filecmp.cmp(new + f, prod + f, shallow=False)}")
            print("  W/U per aperture: " + ' '.join(f'{x:.3f}' for x in ts_w / ts_u))


if __name__ == '__main__':
    main()
