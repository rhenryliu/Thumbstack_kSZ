"""Whole-range against per-sub-bin velocity zeroing for the BGS_BRIGHT-21.35 bin 4 (0.1 < z <= 0.4) kSZ stacks of
configs/velocity_dr2/BGS_BRIGHT-21.35_dsigma_lensing_*.yaml (read-only).

Recomputes ThumbStack's ksz_uniformweight estimator (PS mask and 5-sigma-PTE outlier rejection as catalogMask; saved
profiles use the catalogue's r_V = 1) on the shared bin-4 filtering with (a) the bin-4 velocities, zeroed per cap over
the whole range as stacked, and (b) the per-cap, per-sub-bin zeroed velocities of bins 1-3 for the same objects.
(a) must reproduce the saved profiles. Prints the residual per-cap, per-sub-range mean velocity in bin 4 and the simple
template ratio sum(x o) / sum(o o) of each set to the official one for (a) and (b).

Recorded run: 2026-11-19, on the stacks of interactive jobs 61608703 and 61611035: (a) reproduces the saved profiles to
2e-14; residual means up to +-113 km/s; ratios to official, (a) -> (b): boryana_MG 1.066 -> 0.909, v_recsym
1.061 -> 1.022, v_lift_R5 1.024 -> 0.859; the official profile itself moves by 30% at R = 1 arcmin.

Usage: python bgs_bin4_zeroing.py   (login node, cosmodesi environment)
"""
import numpy as np
import pandas as pd
from scipy import optimize, special

CAT = '/pscratch/sd/r/rhliu/projects/Weak_lensing/desi/spec_Y3/BGS_catalogues_zeroV_vel_dr2/BGS_BRIGHT-21.35/'
OUT = '/pscratch/sd/r/rhliu/projects/Weak_lensing/ksz_measurements/ACTxDESI/spec_Y3_DSigma_lensing_vel_dr2/BGS_BRIGHT-21.35/'
F = OUT + 'official/full_no_src_with_cluster_mask_zeroV/DESIY3_BGS_z4/'
t = np.loadtxt(F + 'DSigma_lensing_filtmap.txt')
fm = np.loadtxt(F + 'DSigma_lensing_filtmask.txt')
mask = np.abs(fm[:, -1]) < 1.
# ThumbStack's outlier rejection (catalogMask): drop objects beyond the 5-sigma-PTE threshold in any aperture
n = mask.sum()
cut = optimize.brentq(lambda x: n * special.erfc(x / np.sqrt(2.)) - special.erfc(5. / np.sqrt(2.)), 0., 1.e2)
mask &= np.all(np.abs(t) <= cut * np.std(t[mask], axis=0)[None, :], axis=1)
print(f"bin 4: {len(t)} objects, {mask.sum()} kept by the PS mask and outlier rejection")
cols = ['TARGETID', 'Z', 'GC', 'VEL_LOS_RENORM']


def stack(vel):
    v = -(vel - vel.mean())[mask] / 3e5 # stack-time zeroV over the catalogue, then -v/c as ThumbStack
    return np.std(v) * (v @ t[mask]) / (v @ v)


prof = {}
for s in ['official', 'boryana_MG', 'v_recsym', 'v_lift_R5']:
    d4 = pd.read_csv(CAT + s + '/full_zbin4_no_src_with_cluster_mask.csv', usecols=cols, float_precision='round_trip')
    sub = pd.concat([pd.read_csv(CAT + s + f'/full_zbin{b}_no_src_with_cluster_mask.csv', usecols=cols,
                                 float_precision='round_trip') for b in (1, 2, 3)]).set_index('TARGETID')
    v_a = d4['VEL_LOS_RENORM'].values
    v_b = sub['VEL_LOS_RENORM'].reindex(d4['TARGETID'].values).values
    assert np.isfinite(v_b).all()
    meas = np.loadtxt(OUT + s + '/full_no_src_with_cluster_mask_zeroV/DESIY3_BGS_z4/DSigma_lensing_ksz_uniformweight_measured.txt')[:, 1]
    pa, pb = stack(v_a), stack(v_b)
    prof[s] = (pa, pb)
    print(f"\n{s}: reproduce measured (saved profiles use r_V = 1): max |pa / meas - 1| = {np.max(np.abs(pa / meas - 1)):.2e}")
    z = d4['Z'].values
    for lo, hi in ((0.1, 0.2), (0.2, 0.3), (0.3, 0.4)):
        for g in ('NGC', 'SGC'):
            sel = (z > lo) & (z <= hi) & (d4['GC'].values == g)
            print(f"   {lo}-{hi} {g}: residual mean of bin-4 velocity {v_a[sel].mean():+7.1f} km/s (std {v_a[sel].std():.0f})")
for s in ['boryana_MG', 'v_recsym', 'v_lift_R5']:
    for k, lab in ((0, 'a: bin-4 zeroing (as stacked)'), (1, 'b: per-sub-bin zeroing')):
        x, o = prof[s][k], prof['official'][k]
        print(f"{s}/official {lab}: simple amplitude sum(x o)/sum(o o) = {(x @ o) / (o @ o):.3f}")
for k, lab in ((0, 'a'), (1, 'b')):
    print(f"profiles ({lab}) official: {np.round(prof['official'][k], 4)}  lift: {np.round(prof['v_lift_R5'][k], 4)}")
