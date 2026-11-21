"""Floor masses M_min of the stellar-mass weights w = M*^p + M_min^p (configs/mstarweight_grid/), one per sample.

M_min is the mean stellar mass (linear, not log) of the q lowest-mass fraction of the full-sky catalogue,
objects with LOGMSTAR <= 0 (failed fits) excluded; printed as log10 M_min. The NGC and SGC stacks of LRG z2
use the full-sky z2 value. The values go into the configs to 4 decimals.

Usage: python mstarweight_floor_values.py  (login node; reads one column of each catalogue)
Recorded run: 2026-11-20, the values in configs/mstarweight_grid/.
"""
import numpy as np
import pandas as pd

W = '/pscratch/sd/r/rhliu/projects/Weak_lensing/desi/spec_Y3/'
CUT = 'no_src_with_cluster_mask'
CATS = {f'DESIY3_LRG_z{b}': W + f'LRG_catalogues_zeroV/full_zbin{b}_{CUT}.csv' for b in (1, 2, 3, 4)}
CATS['DESIY3_BGS'] = W + f'BGS_CIGALE/BGS_catalogues/full_catalog_Y3_{CUT}.csv'
QUANTILES = (0.05, 0.10, 0.25)


def log_mean_lowest(logm: np.ndarray, q: float) -> float:
    """log10 of the mean of 10^logm over the objects at or below the q quantile of logm."""
    low = logm[logm <= np.quantile(logm, q)]
    return float(np.log10(np.mean(10.**low)))


def main():
    for q in QUANTILES:
        print(f"floor{round(100 * q):02d}: log10 M_min, mean M* of the lowest {100 * q:.0f} %")
        for name, path in CATS.items():
            logm = pd.read_csv(path, usecols=['LOGMSTAR'])['LOGMSTAR'].values.astype(float)
            logm = logm[np.isfinite(logm) & (logm > 0.)]
            print(f"  {name}: {log_mean_lowest(logm, q):.4f}  ({len(logm)} objects)")


if __name__ == '__main__':
    main()
