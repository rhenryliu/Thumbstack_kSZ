import sys
sys.path.append('/global/homes/r/rhliu/projects/repos/ThumbStack/')

from importlib import reload
import universe
reload(universe)
from universe import *

import mass_conversion
reload(mass_conversion)
from mass_conversion import *

import catalog
reload(catalog)
from catalog import *

import thumbstack
reload(thumbstack)
from thumbstack import *

import cmb
reload(cmb)
from cmb import *
from cmbMap import *
import matplotlib as mpl
mpl.rcParams.update(mpl.rcParamsDefault)
import json

from catalog_utils import make_Catalog
import pandas as pd

# This script follows prepare_catalog.py to convert DESI DR2 catalogues to ThumbStack format.
# The main reason this script is separate from prepare_catalogue is that prepare_catalogue (written by Frank)
# requires a different python environment to run. This one requires rhliu_tSZ (my specific ThumbStack environment).

pathCat = '/pscratch/sd/r/rhliu/projects/Weak_lensing/desi/spec_Y3/'
data_names = ['RA', 'DEC', 'Z', 'VEL_LOS_RENORM']

# data = np.loadtxt(pathCat + 'full_zbin1_unfiltered.txt')
# cat = {name: data[:, i] for i, name in enumerate(data_names)}

bins = [1, 2, 3, 4]

u = UnivMariana()
massConversion = MassConversionKravtsov14()
out_dir = '/pscratch/sd/r/rhliu/projects/Weak_lensing/desi/spec_Y3/'
fig_dir = '/pscratch/sd/r/rhliu/projects/Weak_lensing/figs/'

field = 'full'  # 'NGC' or 'SGC' or 'full'
fields = ['NGC', 'SGC', 'full']
filter_type = 'unfiltered'  # 'unfiltered' or 'no_src_with_cluster_mask'
# filter_type = 'no_src_with_cluster_mask'


# cat_fn = "/catalog.txt"
# cat_fn = "/catalog_no_src_with_cluster_mask.txt"
# cat_fn = "/catalog_NGC_unfiltered.txt"
# cat_fn = "/catalog_SGC_unfiltered.txt"
# cat_fn = f"/catalog_{field}_{filter_type}.txt"

for field in fields:
    print(f'Processing field {field}...')
    cat_fn = f"/catalog_{field}_{filter_type}.txt"
    # data_df = pd.read_csv(pathCat + f'full_zbin{bin}_unfiltered.csv')
    # data_df = pd.read_csv(pathCat + f'full_zbin{bin}_no_src_with_cluster_mask.csv')
    # data_df = pd.read_csv(pathCat + f'NGC_zbin{bin}_unfiltered.csv')
    data_df = pd.read_csv(pathCat + f'{field}_catalog_Y3_{filter_type}.csv')
    data_df = data_df[data_names]
    data_df.rename(columns={'VEL_LOS_RENORM': 'vR'}, inplace=True) # type: ignore
    print(f'Dataframe shape: {data_df.shape}')

    cat_name = f'DESIY3_LRG'
    
    cat = make_Catalog(u, massConversion, data_df, 
                       out_dir=out_dir,
                       fig_dir=fig_dir, 
                       cat_fn=cat_fn,
                       name=cat_name)
    print('Saving to ', out_dir + cat_name + cat_fn)
    cat.writeCatalog()
    
print('All done!')
# data_df = pd.read_csv(pathCat + 'full_zbin1_unfiltered.csv')
# data_df = data_df[data_names]
# data_df.rename(columns={'VEL_LOS_RENORM': 'vR'}, inplace=True) # type: ignore

