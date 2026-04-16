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
# Velocity-shuffle null test with multiple seeds.
#
# Relevant YAML parameters (under 'stack'):
#   shuffle_seed : int   — master seed used to generate the list of shuffle seeds
#   n_seeds      : int   — number of shuffle realizations to run
#   bootstrap    : bool  — if True, each realization runs a full bootstrap to
#                          estimate its covariance; combined covariance is
#                          (1/N) * mean(Cov_i).
#                          if False, covariance is estimated from the inter-shuffle
#                          sample covariance divided by N (no bootstrap at all,
#                          much faster; requires n_seeds >> nRAp for reliability).
#
# Output files are the same format as make_stacks_yaml.py:
#   {filterType}_{est}_measured.txt          — [R, mean_signal, combined_error]
#   cov_{filterType}_{est}_bootstrap.txt     — combined NxN covariance matrix
##################################################################################
t0 = time.time()

extra_str = "" # initialize
want_sim = False

parser = argparse.ArgumentParser(description='Process config.')
parser.add_argument('-p', '--path2config', type=str, default='./configs/null/LRG_dsigma_shuffle_zeroV.yaml', help='Path to the configuration file.')
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
cat_names = args['cat_names'] if args['cat_names'] is not None else data_config.get('cat_names', [f"DESIY3_{cat_type}"]) # Names for the different catalogues, which will be used to read in the correct files and save the correct output files. Should be a list of strings.

cat_dir = data_config.get('cat_dir', '/pscratch/sd/r/rhliu/projects/Weak_lensing/desi/spec_Y3/{cat_type}_catalogue/').format(cat_type=cat_type) # Frank's SPEC Y3 data
fig_dir = data_config.get('fig_dir', '/pscratch/sd/r/rhliu/projects/Weak_lensing/figs/')
cat_fn_template = data_config.get('cat_fn', 'catalog_{cat_type}_{field}_{filter_cut}.txt')
cat_fn = cat_fn_template.format(field=field, filter_cut=filter_cut)
output_dir_ = data_config.get('output_dir', f"/pscratch/sd/r/rhliu/projects/Weak_lensing/ksz_measurements/ACTxDESI/spec_Y3/")
extra_str = data_config.get('extra_str', f"{field}_{filter_cut}").format(field=field, filter_cut=filter_cut)

# --- SHUFFLE NULL TEST PARAMETERS ---
shuffle_seed = stack_config.get('shuffle_seed', 42)  # master seed for generating the seed list
n_seeds = stack_config.get('n_seeds', 10)            # number of shuffle realizations
doBootstrap = stack_config.get('bootstrap', True)    # covariance method (see module docstring above)

# Generate n_seeds shuffle seeds from the master seed
_seed_rng = np.random.default_rng(seed=shuffle_seed)
shuffle_seeds = _seed_rng.integers(0, 2**31, size=n_seeds).tolist()
print(f"Shuffle null test: {n_seeds} realizations, master seed {shuffle_seed}")
print(f"Generated seeds: {shuffle_seeds}")
print(f"Covariance method: {'bootstrap (per-seed)' if doBootstrap else 'inter-shuffle sample covariance'}")

if not doBootstrap and n_seeds < 20:
    print(f"WARNING: inter-shuffle covariance with n_seeds={n_seeds} may be unreliable. "
          f"Recommend n_seeds >= 50 for a stable {stack_config.get('filterType','DSigma')} covariance matrix.")

# cosmological parameters
u = UnivMariana()
massConversion = MassConversionKravtsov14()
invPowerFunc = None
filterFuncRad = None
apod_pix = 20

# Stack config params;
save = stack_config.get('save', False) # if yes calculates delta T anew (saves mask, delta T, maybe vel is separate? but the numbers need to match)
doMBins = False # we don't have mass bins yet
doVShuffle = False #True # default measurement is False
wantMF = False # default measurement is False
doOnlyFiltering = False #True # default is False; if True, compute dT decrements and return TESTING!!!!!!!!!!!; the only problem is that.... it doesn't compute covariance only stackedMap; we use this for anisotropic

mode = "kSZ"
if "tau_screening" in mode:
    filterType = "meanring" # og
    Obs = 'tau'
elif "lensing" in mode:
    filterType = "meanring"
    Obs = 'tsz'
else:
    filterType = filterType # use the filter type from the config, which can be overridden by command line argument
    Obs = 'ksz' # og

if stack_config.get('zeroV', False):
    print("Removing mean velocity from the catalog")
    extra_str += "_zeroV"

pathMap = "/pscratch/sd/b/boryanah/ACTxDESI/ACT/hilc_fullRes_TT_17000.fits" # 1.6 arcmin # OG
pathMask = '/pscratch/sd/b/boryanah/ACTxDESI/ACT/wide_mask_GAL070_apod_1.50_deg_wExtended_srcfree_Will.fits'
output_dir = output_dir_ + f"{extra_str}/"
pathHit = None

catalogs = {}
for i, cat_name in enumerate(cat_names):
    catalogs[cat_name] = Catalog(u, massConversion, name=cat_name, nameLong=cat_name, out_dir=cat_dir, save=False, fig_dir=fig_dir, cat_fn=cat_fn)

# Read CMB maps
nProc = 128

CMB_nu = 90.e9
cmap = cmbMap(pathMap,
              pathMask=pathMask,
              pathHit=pathHit,
              nu=CMB_nu, unitLatex=r'y',
              name='')

# --- MASS FILTERING SETUP (optional) ---
# Read once here; per-catalog application happens inside the loop below.
# If 'mass' is absent from the config this is a no-op and all existing
# behaviour is preserved.
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
        _mass_mask = compute_mass_mask(_logmstar, _mass_config) # type: ignore
        _n_before = catalog.nObj
        apply_catalog_mask(catalog, _mass_mask)
        print(f"[mass filter] {key}: {_n_before} -> {catalog.nObj} objects retained")

    # Save original velocities; restored before each seed so every shuffle
    # starts from the same unmodified vR.
    vR_orig = catalog.vR.copy()

    # Accumulators across shuffle seeds
    all_signals = {}  # {filterType_est: list of arrays shape (nRAp,)}
    all_covs = {}     # {filterType_est: list of arrays shape (nRAp, nRAp)}, only when doBootstrap=True

    for i_seed, seed in enumerate(shuffle_seeds):
        print(f"\n--- Catalog {key}: shuffle realization {i_seed+1}/{n_seeds} (seed={seed}) ---")

        # Restore original velocities for this seed
        catalog.vR = vR_orig.copy()

        #### CHANGE HERE IF NECESSARY: ZEROING VELOCITY COMPONENT:
        if stack_config.get('zeroV', False):
            print("Removing mean velocity from the catalog")
            catalog.vR -= np.mean(catalog.vR) # removing mean velocity # type: ignore

        print(f"Shuffling velocities in the catalog for null test (seed={seed})")
        rng = np.random.default_rng(seed=seed)
        rng.shuffle(catalog.vR) # shuffling velocity for null test # type: ignore

        if Obs == 'tau':
            catalog.vR = catalog.vZ # hiding here info about T_large-scales

        # Compute and save the CMB temperature decrements (filtMap) only on the
        # first seed; all subsequent seeds load the already-saved filtMap.
        # This avoids redundant cutout extraction, saving significant compute time.
        save_this_seed = save if i_seed == 0 else False

        ts = ThumbStack(u, catalog,
                        cmap.map(),
                        cmap.mask(),
                        cmap.hit(),
                        name=catalog.name,
                        nameLong=catalog.nameLong,
                        save=save_this_seed,
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
                        apod_pix=apod_pix)

        # Collect results from this seed's ThumbStack.
        # ts.stackedProfile and ts.covBootstrap are populated by
        # loadAllStackedProfiles() called inside the ThumbStack constructor.
        for ft in ts.filterTypes:
            for est in ts.Est:
                result_key = ft + "_" + est
                if result_key not in all_signals:
                    all_signals[result_key] = []
                    all_covs[result_key] = []
                all_signals[result_key].append(ts.stackedProfile[result_key].copy())
                if doBootstrap:
                    all_covs[result_key].append(ts.covBootstrap[result_key].copy())

    # --- Combine results across all seeds and write output files ---
    print(f"\n--- Catalog {key}: combining {n_seeds} shuffle realizations ---")
    for ft in ts.filterTypes:
        for est in ts.Est:
            result_key = ft + "_" + est
            signals_arr = np.array(all_signals[result_key])  # (n_seeds, nRAp)
            mean_signal = np.mean(signals_arr, axis=0)       # (nRAp,)

            if doBootstrap:
                # Covariance of the mean of N independent measurements:
                # Cov(mean) = (1/N^2) * sum(Cov_i) = (1/N) * mean(Cov_i)
                covs_arr = np.array(all_covs[result_key])    # (n_seeds, nRAp, nRAp)
                combined_cov = np.mean(covs_arr, axis=0) / n_seeds
            else:
                # Inter-shuffle sample covariance: estimate variance of the null
                # test statistic from the scatter across shuffle realizations.
                # Covariance of the mean is the sample covariance divided by N.
                combined_cov = np.cov(signals_arr, rowvar=False) / n_seeds

            combined_error = np.sqrt(np.diag(combined_cov))

            # Overwrite the on-disk files with the combined result
            data = np.zeros((ts.nRAp, 3))
            data[:, 0] = ts.RApArcmin
            data[:, 1] = mean_signal
            data[:, 2] = combined_error
            np.savetxt(ts.pathOut + "/" + ft + "_" + est + "_measured.txt", data)
            print(f"Saved combined measured profile: {ts.pathOut}/{ft}_{est}_measured.txt")

            np.savetxt(ts.pathOut + "/cov_" + ft + "_" + est + "_bootstrap.txt", combined_cov)
            print(f"Saved combined covariance:       {ts.pathOut}/cov_{ft}_{est}_bootstrap.txt")

    # Now we want to add some code to save the individual stacks.
    # Note: for the multi-seed shuffle, this block uses the last seed's ThumbStack (ts).
    if extra_config.get('save_individual_stacks', False):
    # if False: # for now we don't need to save individual stacks, and this is a bit slow, so I'm commenting it out. We can uncomment and use when needed.

        # Note that the catalogues here are the original catalogues before the processing necessary to read into ThumbStack
        # (e.g. adding columns), so we read in the original catalogues again to get the original vR values (before zeroV correction, if applied).
        # We also add the stacks as new columns to this original catalogue and save it as a new csv file.

        catalogue_path = extra_config.get('catalogue_path', '/pscratch/sd/r/rhliu/projects/Weak_lensing/desi/spec_Y3/{cat_type}_catalogues/').format(cat_type=cat_type)
        # catalogue_dataframe_template = extra_config.get('catalogue_dataframe', '{field}_{bin}_{filter_cut}.csv')
        catalogue_dataframe_template = '{field}_{bin}_{filter_cut}.csv'

        if key in ["DESIY3_LRG", "DESIY3_BGS"]:
            # for the all z bins catalogues, we use a different naming convention for the dataframe, since they don't correspond to a specific z bin.
            catalogue_dataframe = catalogue_dataframe_template.format(field=field, bin='catalog_Y3', filter_cut=filter_cut)
        else:
            sample, zbin = key.rsplit('_z', 1) # split on the last occurrence of '_z' to separate the sample name and the z bin number
            catalogue_dataframe = catalogue_dataframe_template.format(field=field, bin=f'zbin{zbin}', filter_cut=filter_cut)

        # df_catalog = pd.read_csv(cat_dir + cat_fn, delim_whitespace=True)
        data_names = ['TARGETID', 'Z','RA', 'DEC', 'VEL_LOS_RENORM']
        data_names = ['Z','RA', 'DEC', 'VEL_LOS_RENORM']
        data_df = pd.read_csv(catalogue_path + catalogue_dataframe)
        data_df = data_df[data_names]
        data_df.rename(columns={'VEL_LOS_RENORM': 'vR'}, inplace=True)

        allProfiles = computeProfiles(ts, filterType, est='ksz_uniformweight')
        print("Computed all profiles for individual stacks. Allprofiles shape:", allProfiles.shape)
        mask = ts.catalogMask(overlap=True, psMask=True, filterType=filterType, mVir=None)

        # data_df['Mask'] = mask.astype(int)
        data_df = data_df[mask]
        NObj = data_df.shape[0]

        for j, R in enumerate(ts.RApArcmin):
            print(str(R), )
            # stacks = np.zeros(NObj)
            # stacks[mask] = allProfiles[:, j]
            stacks = allProfiles[:, j]

            data_df[str(R)] = stacks

        # save_name = extra_config.get('save_dir', '').format(field=field, filter_cut=filter_cut) + extra_config.get('save_name', '{filterType}_profiles_{field}_{filter_cut}_z{bin}.csv').format(filterType=filterType, field=field, bin=i+1, filter_cut=filter_cut)
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
