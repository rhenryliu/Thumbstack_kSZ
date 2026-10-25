"""Known answer for the new DSigma_lensing cap filtering (one-off check, kept as the record behind docs/completeness_weighted_ksz_results.md): the per-object
filtering of each NGC/SGC object must equal that object's row in the production full-sky files
(same map, same position), objects matched by their exact RA and DEC.

Usage: python cap_filtering_check.py LRG_z1 [LRG_z4 BGS ...]
Recorded runs: LRG_z1 to LRG_z3 (cap_filtering_check_LRG.txt), BGS (job 60185347 log).
"""
import sys

import numpy as np
import pandas as pd

D = '/pscratch/sd/r/rhliu/projects/Weak_lensing/desi/spec_Y3/'
K = '/pscratch/sd/r/rhliu/projects/Weak_lensing/ksz_measurements/ACTxDESI/'
CUT = 'no_src_with_cluster_mask'
FILES = ('filtmap', 'filtmask', 'filtnoisestddev', 'filtarea')
SAMPLES = {f'LRG_z{b}': dict(cat=D + f'LRG_catalogues_zeroV/DESIY3_LRG_z{b}/catalog_{{f}}_{CUT}.txt',
                             full_cat=D + f'LRG_catalogues_zeroV/DESIY3_LRG_z{b}/catalog_full_{CUT}.txt',
                             full=K + f'spec_Y3_DSigma_lensing/full_{CUT}_zeroV/DESIY3_LRG_z{b}/',
                             cap=K + f'spec_Y3_DSigma_lensing_compweight/{{f}}_{CUT}_zeroV/DESIY3_LRG_z{b}/') for b in (1, 2, 3, 4)}
# caps from the non-CIGALE catalogue; the full-sky v2-filter production stack is the CIGALE one
SAMPLES['BGS'] = dict(cat=D + f'BGS_catalogues/DESIY3_BGS/catalog_{{f}}_{CUT}.txt',
                      full_cat=D + f'BGS_CIGALE/BGS_catalogues/DESIY3_BGS/catalog_full_{CUT}.txt',
                      full=K + f'spec_Y3_DSigma_lensing/CIGALE/full_{CUT}_zeroV/DESIY3_BGS/',
                      cap=K + f'spec_Y3_DSigma_lensing_compweight/{{f}}_{CUT}_zeroV/DESIY3_BGS/')


def radec(path):
    return pd.read_csv(path, sep=r'\s+', header=None, usecols=[0, 1]).values


def table(path):
    return pd.read_csv(path, sep=r'\s+', header=None).values


for name in sys.argv[1:]:
    s = SAMPLES[name]
    full_pos = radec(s['full_cat'])
    index = {tuple(p): i for i, p in enumerate(full_pos)}
    full = {f: table(s['full'] + f'DSigma_lensing_{f}.txt') for f in FILES}
    for field in ('NGC', 'SGC'):
        pos = radec(s['cat'].format(f=field))
        rows = np.array([index.get(tuple(p), -1) for p in pos])
        found = rows >= 0
        diffs = []
        for f in FILES:
            cap = table(s['cap'].format(f=field) + f'DSigma_lensing_{f}.txt')
            diffs.append(np.max(np.abs(cap[found] - full[f][rows[found]])))
        print(f"{name} {field}: {found.sum()} of {len(pos)} objects found in the full-sky catalogue; "
              f"max |cap - full| over {', '.join(FILES)}: {', '.join(f'{d:.1e}' for d in diffs)}")
