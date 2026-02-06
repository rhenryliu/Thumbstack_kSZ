import sys

import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from pathlib import Path
from datetime import datetime

# This script follows prepare_catalog.py to convert DESI DR2 catalogues to ThumbStack format.
# The main reason this script is separate from prepare_catalogue is that prepare_catalogue (written by Frank)
# requires a different python environment to run. This one requires rhliu_tSZ (my specific ThumbStack environment).

pathCat = '/pscratch/sd/r/rhliu/projects/Weak_lensing/desi/spec_Y3/'
field = 'NGC'
filter_ = 'unfiltered'
# catfile = f'{field}_zbin{{bin}}_no_src_with_cluster_mask.csv'
catfile = f'{field}_zbin{{bin}}_{filter_}.csv'
data_names = ['RA', 'DEC', 'Z', 'VEL_LOS_RENORM']

now = datetime.now()
yr_string = now.strftime("%Y-%m")
dt_string = now.strftime("%m-%d")
figPath = Path('../figures/') / yr_string / dt_string
figPath.mkdir(parents=True, exist_ok=True)

# data = np.loadtxt(pathCat + 'full_zbin1_unfiltered.txt')
# cat = {name: data[:, i] for i, name in enumerate(data_names)}

bins = [1, 2, 3, 4]

# out_dir = '/pscratch/sd/r/rhliu/projects/Weak_lensing/desi/spec_Y3/'
# fig_dir = '/pscratch/sd/r/rhliu/projects/Weak_lensing/figs/'
# cat_fn = "/catalog.txt"
# cat_fn = "/catalog_no_src_with_cluster_mask.txt"
# cat_fn = "/catalog_NGC_unfiltered.txt"
# cat_fn = "/catalog_SGC_unfiltered.txt"
# cat_fn = "/catalog_NGC_no_src_with_cluster_mask.txt"


fig, axes = plt.subplots(2, 2, figsize=(10, 8))
axes = axes.flatten()
for i, bin in enumerate(bins):
    ax = axes[i]
    print(f'Processing z bin {bin}...')
    # data_df = pd.read_csv(pathCat + f'full_zbin{bin}_unfiltered.csv')
    # data_df = pd.read_csv(pathCat + f'full_zbin{bin}_no_src_with_cluster_mask.csv')
    # data_df = pd.read_csv(pathCat + f'NGC_zbin{bin}_unfiltered.csv')
    data_df = pd.read_csv(pathCat + catfile.format(bin=bin))
    # data_df = data_df[data_names]
    
    
    # data_df.rename(columns={'VEL_LOS_RENORM': 'vR'}, inplace=True) # type: ignore
    print(f'Dataframe shape: {data_df.shape}')

    cat_name = f'DESIY3_LRG_z{bin}, N={data_df.shape[0]}, mean vel={np.mean(data_df["VEL_LOS_RENORM"]):.2f} km/s'
    
    ax.hist(data_df['VEL_LOS_RENORM'], bins=50, alpha=0.7, density=True)
    ax.axvline(np.mean(data_df['VEL_LOS_RENORM']), color='r', linestyle='dashed', linewidth=1)
    ax.set_title(cat_name)
    ax.set_xlabel('VEL_LOS_RENORM')
    ax.set_ylabel('Number of galaxies')
    ax.grid()
plt.tight_layout()
plt.savefig(figPath / f'desi_specY3_velocities_{field}_{filter_}.png')
    
    
print('All done!')
# data_df = pd.read_csv(pathCat + 'full_zbin1_unfiltered.csv')
# data_df = data_df[data_names]
# data_df.rename(columns={'VEL_LOS_RENORM': 'vR'}, inplace=True) # type: ignore

