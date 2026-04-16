"""Median-estimator diagnostic for kSZ stacking.

Identical to make_stacks_yaml.py except:
  - Reads `stack.median` (default True) and passes use_median to ThumbStack.
  - Reads `stack.bootstrap` (default False) for doBootstrap.
  - Appends `_median` to extra_str when median=True, routing output to a
    separate directory so existing mean-estimator results are never overwritten.

Run exactly like the main script:
    python -u make_stacks_median.py -p ./configs/null/LRG_dsigma_median.yaml
    python -u make_stacks_median.py -p ./configs/null/LRG_dsigma_median.yaml \
        --filterType DSigma --field NGC --filter-cut unfiltered
"""
import sys
sys.path.append('/global/homes/r/rhliu/projects/repos/ThumbStack_kSZ/')

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
import yaml
import argparse
import pandas as pd
from computeProfiles import computeProfiles
import time
import os

##################################################################################
t0 = time.time()

extra_str = "" # initialize
want_sim = False

parser = argparse.ArgumentParser(description='Process config.')
parser.add_argument('-p', '--path2config', type=str, default='./configs/null/LRG_dsigma_median.yaml', help='Path to the configuration file.')
parser.add_argument('--filterType', type=str, default=None, help='Override filterType from config (DSigma or diskring).')
parser.add_argument('--field', type=str, default=None, help='Override field from config (NGC, SGC, or full).')
parser.add_argument('--filter-cut', type=str, default=None, help='Override filter_cut from config.')
parser.add_argument('--cat-type', type=str, default=None, help='Override cat_type from config (LRG or BGS).')
parser.add_argument('--cat-names', type=str, nargs='+', default=None, help='Override cat_names from config (space-separated list).')
args = vars(parser.parse_args())
print(f"Arguments: {args}")
path2config = args['path2config']
with open(path2config, 'r') as f:
    config = yaml.safe_load(f)

stack_config = config.get('stack', {})
data_config = config.get('data', {})
extra_config = config.get('extra', {})

filterType = args['filterType'] if args['filterType'] is not None else stack_config.get("filterType", "DSigma")
cat_type = args['cat_type'] if args['cat_type'] is not None else data_config.get('cat_type', 'LRG')
field = args['field'] if args['field'] is not None else data_config.get('field', 'NGC')
filter_cut = args['filter_cut'] if args['filter_cut'] is not None else data_config.get('filter_cut', 'unfiltered')
cat_names = args['cat_names'] if args['cat_names'] is not None else data_config.get('cat_names', [f"DESIY3_{cat_type}"])

cat_dir = data_config.get('cat_dir', '/pscratch/sd/r/rhliu/projects/Weak_lensing/desi/spec_Y3/{cat_type}_catalogue/').format(cat_type=cat_type)
fig_dir = data_config.get('fig_dir', '/pscratch/sd/r/rhliu/projects/Weak_lensing/figs/')
cat_fn_template = data_config.get('cat_fn', 'catalog_{cat_type}_{field}_{filter_cut}.txt')
cat_fn = cat_fn_template.format(field=field, filter_cut=filter_cut)
output_dir_ = data_config.get('output_dir', f"/pscratch/sd/r/rhliu/projects/Weak_lensing/ksz_measurements/ACTxDESI/spec_Y3/")
extra_str = data_config.get('extra_str', f"{field}_{filter_cut}").format(field=field, filter_cut=filter_cut)

# cosmological parameters
u = UnivMariana()
massConversion = MassConversionKravtsov14()
invPowerFunc = None
filterFuncRad = None
apod_pix = 20

# Stack config params
save = stack_config.get('save', False)
doMBins = False
doVShuffle = False
wantMF = False
doOnlyFiltering = False

# Median diagnostic flag: use median instead of mean in the stacked profile estimator.
# Defaults to True since this script is specifically for the median diagnostic.
use_median = stack_config.get('median', True)

# Bootstrap is expensive and the approximate sStack is sufficient for this diagnostic.
doBootstrap = stack_config.get('bootstrap', False)

mode = "kSZ"
if "tau_screening" in mode:
    filterType = "meanring"
    Obs = 'tau'
elif "lensing" in mode:
    filterType = "meanring"
    Obs = 'tsz'
else:
    filterType = filterType
    Obs = 'ksz'

if stack_config.get('zeroV', False):
    print("Removing mean velocity from the catalog")
    extra_str += "_zeroV"

# Append _median suffix so outputs land in a separate directory,
# never overwriting existing mean-estimator results.
if use_median:
    extra_str += "_median"
    print(f"Median estimator enabled: output will be written to .../{extra_str}/")

pathMap = "/pscratch/sd/b/boryanah/ACTxDESI/ACT/hilc_fullRes_TT_17000.fits"
pathMask = '/pscratch/sd/b/boryanah/ACTxDESI/ACT/wide_mask_GAL070_apod_1.50_deg_wExtended_srcfree_Will.fits'
output_dir = output_dir_ + f"{extra_str}/"
pathHit = None

catalogs = {}
for i, cat_name in enumerate(cat_names):
    catalogs[cat_name] = Catalog(u, massConversion, name=cat_name, nameLong=cat_name, out_dir=cat_dir, save=False, fig_dir=fig_dir, cat_fn=cat_fn)

nProc = 128

CMB_nu = 90.e9
cmap = cmbMap(pathMap,
              pathMask=pathMask,
              pathHit=pathHit,
              nu=CMB_nu, unitLatex=r'y',
              name='')

# --- MASS FILTERING SETUP (optional) ---
_mass_config = config.get('mass', None)
if _mass_config is not None:
    from mass_filter_utils import (compute_mass_mask, apply_catalog_mask,
                                   get_csv_path, mass_label as _mass_label_fn)
    _mlabel = _mass_label_fn(_mass_config)
    effective_output_dir = output_dir_ + f"{extra_str}_{_mlabel}/"
    print(f"[mass filter] Active: strategy='{_mass_config['strategy']}', "
          f"label='{_mlabel}'")
    print(f"[mass filter] Output directory: {effective_output_dir}")
else:
    _mlabel = None
    effective_output_dir = output_dir

for i, key in enumerate(catalogs.keys()):
    catalog = catalogs[key]
    catalog.Mstellar = np.empty_like(catalog.RA)
    catalog.Mvir = np.empty_like(catalog.RA)
    catalog.integratedY = np.empty_like(catalog.RA)
    catalog.integratedKSZ = np.empty_like(catalog.RA)
    catalog.integratedTau = np.empty_like(catalog.RA)

    # --- MASS FILTERING (per catalog) ---
    if _mass_config is not None:
        _mass_col = _mass_config.get('col', 'LOGMSTAR')
        _csv_path = get_csv_path(key, cat_type, cat_dir, field, filter_cut)
        print(f"[mass filter] {key}: reading '{_mass_col}' from {_csv_path}")
        _mass_df = pd.read_csv(_csv_path, usecols=[_mass_col])
        _logmstar = _mass_df[_mass_col].values
        if len(_logmstar) != catalog.nObj:
            raise ValueError(
                f"[mass filter] Mass column length ({len(_logmstar)}) != "
                f"catalog nObj ({catalog.nObj}) for {key}. "
                f"CSV: {_csv_path}. Ensure the CSV and ThumbStack catalog "
                f"were produced from the same run."
            )
        _mass_mask = compute_mass_mask(_logmstar, _mass_config)
        _n_before = catalog.nObj
        apply_catalog_mask(catalog, _mass_mask)
        print(f"[mass filter] {key}: {_n_before} -> {catalog.nObj} objects retained")

    if stack_config.get('zeroV', False):
        print("Removing mean velocity from the catalog")
        catalog.vR -= np.mean(catalog.vR)

    if stack_config.get('shuffle', False):
        print("Shuffling velocities in the catalog for null test")
        seed = stack_config.get('shuffle_seed', 42)
        rng = np.random.default_rng(seed=seed)
        rng.shuffle(catalog.vR)

    if Obs == 'tau':
        catalog.vR = catalog.vZ

    ts = ThumbStack(u, catalog,
                    cmap.map(),
                    cmap.mask(),
                    cmap.hit(),
                    name=catalog.name,
                    nameLong=catalog.nameLong,
                    save=save,
                    nProc=nProc,
                    filterTypes=filterType,
                    doMBins=doMBins,
                    doBootstrap=doBootstrap,
                    doVShuffle=doVShuffle,
                    doOnlyFiltering=doOnlyFiltering,
                    cmbNu=cmap.nu,
                    cmbUnitLatex=cmap.unitLatex,
                    output_dir=effective_output_dir,
                    Obs=Obs,
                    wantMF=wantMF,
                    invPowerFunc=invPowerFunc,
                    filterFuncRad=filterFuncRad,
                    apod_pix=apod_pix,
                    use_median=use_median)

    if extra_config.get('save_individual_stacks', False):

        catalogue_path = extra_config.get('catalogue_path', '/pscratch/sd/r/rhliu/projects/Weak_lensing/desi/spec_Y3/{cat_type}_catalogues/').format(cat_type=cat_type)
        catalogue_dataframe_template = '{field}_{bin}_{filter_cut}.csv'

        if key in ["DESIY3_LRG", "DESIY3_BGS"]:
            catalogue_dataframe = catalogue_dataframe_template.format(field=field, bin='catalog_Y3', filter_cut=filter_cut)
        else:
            sample, zbin = key.rsplit('_z', 1)
            catalogue_dataframe = catalogue_dataframe_template.format(field=field, bin=f'zbin{zbin}', filter_cut=filter_cut)

        data_names = ['Z','RA', 'DEC', 'VEL_LOS_RENORM']
        data_df = pd.read_csv(catalogue_path + catalogue_dataframe)
        data_df = data_df[data_names]
        data_df.rename(columns={'VEL_LOS_RENORM': 'vR'}, inplace=True)

        allProfiles = computeProfiles(ts, filterType, est='ksz_uniformweight')
        print("Computed all profiles for individual stacks. Allprofiles shape:", allProfiles.shape)
        mask = ts.catalogMask(overlap=True, psMask=True, filterType=filterType, mVir=None)

        data_df = data_df[mask]
        NObj = data_df.shape[0]

        for j, R in enumerate(ts.RApArcmin):
            print(str(R), )
            stacks = allProfiles[:, j]
            data_df[str(R)] = stacks

        save_dir = extra_config.get('save_dir', '').format(field=field, filter_cut=filter_cut)
        if _mlabel is not None:
            save_dir = save_dir.rstrip('/') + f'_{_mlabel}/'
        save_name = extra_config.get('save_name', '{cat_name}_{filterType}_{field}_{filter_cut}.csv').format(cat_name=key, filterType=filterType, field=field, filter_cut=filter_cut)
        if not os.path.exists(save_dir):
            os.makedirs(save_dir)
        print(f"Saving individual stacks to {save_dir + save_name}")

        data_df.to_csv(save_dir + save_name, index=False)

t1 = time.time()
print(f"Total time: {t1 - t0} seconds")
