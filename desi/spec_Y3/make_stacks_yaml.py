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
# param_Dict = json.loads(sys.argv[1])

##################################################################################
t0 = time.time()

extra_str = "" # initialize
want_sim = False

parser = argparse.ArgumentParser(description='Process config.')
parser.add_argument('-p', '--path2config', type=str, default='./configs/BGS_dsigma_full_cigale.yaml', help='Path to the configuration file.')
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
# extra_str += config.get('extra_str', f"/{cat_type}_{field}_{filter_cut}").format(cat_type=cat_type, field=field, filter_cut=filter_cut)

# if want_sim:
#     #cat_dir = "/global/cfs/cdirs/desi/users/boryanah/kSZ_recon/for_fiona/"
#     #cat_fn = "Extended_LRG_zerr0.0_AbacusSummit_huge_c000_ph201_masked.fits"
#     #cat_fn = "Extended_LRG_zerr2.0_AbacusSummit_huge_c000_ph201_masked.fits"; extra_str += "_zerr2.0"
#     #cat_fn = "Extended_LRG_zerr0.0_AbacusSummit_huge_c000_ph201_masked_tiny.fits"
#     cat_dir = "/pscratch/sd/b/boryanah/websky/abacus/"
#     cat_fn = "halos_LOGMhalo_Msunh13_masked.fits"
# else:
    #cat_dir = "/global/cfs/cdirs/desi/users/boryanah/reconstruction_DESI/recon/"
    # cat_dir = '/pscratch/sd/b/boryanah/ACTxDESI/DESI/' # Boryana's SPEC Y1 data
    # cat_dir = '/pscratch/sd/r/rhliu/projects/Weak_lensing/desi/spec_Y3/' # Frank's SPEC Y3 data
    # fig_dir = '/pscratch/sd/r/rhliu/projects/Weak_lensing/figs/'
    # cat_dir = '/pscratch/sd/b/boryanah/kSZ_pairwise/'
    # cat_dir = '/global/cfs/cdirs/desi/users/boryanah/kSZ_recon/for_fiona/velocities_dr9_extended_lrg_pzbins-dr10_pz_dr10_extended_randoms-1-0-2_remov_isle_nobs1_ebv0.15_stardens2500_lrg_mask_sigmaz0.0500_R12.50_nmesh512_recsym_MG{_bin_1,_bin_2,_bin_3,_bin_4,}_{s,n}gc.npz'
    #perc = 30
    #perc = 10
    #cat_fn = "catalog_BGS_BRIGHT-20.2_R12.50_nmesh512_recsym_MG_masked.fits"
    #cat_fn = f"DESY6_ABSM{perc:d}.fits"
    #cat_fn = "ELG_LOPnotqso_masked.fits"
    #cat_fn = "ELG_LOPnotqso_cigale_masked.fits"
    # cat_fn = "LRG_cigale_masked.fits"
    # cat_fn = "LRG_Y1_cigale_masked.fits"
    # cat_fn = 'LRG.txt'
    # cat_fn = 'catalog_dr10_allfoot_perbin_sigmaz0.0500.txt'

    # cat_fn = 'catalog.txt'
    # These two lines go together, for the filtered catalogues, otherwise use the one above.
    # cat_fn = 'catalog_no_src_with_cluster_mask.txt'
    # extra_str += f"/no_src_with_cluster_mask"
    # 
    # cat_fn = 'catalog_SGC_no_src_with_cluster_mask.txt'
    # extra_str += f"/SGC_no_src_with_cluster_mask"

# cosmological parameters
u = UnivMariana()
massConversion = MassConversionKravtsov14()
invPowerFunc = None
filterFuncRad = None
apod_pix = 20

# Stack config params; 
# TODO: move these to the yaml file if we want to change them more easily
save = stack_config.get('save', False) # if yes calculates delta T anew (saves mask, delta T, maybe vel is separate? but the numbers need to match)
#save = False # False means use saved stuff
doBootstrap = True #False #True # do bootstrap or not? (not needed if e.g. shuffling velocities)
doMBins = False # we don't have mass bins yet
doVShuffle = False #True # default measurement is False
wantMF = False # default measurement is False
doOnlyFiltering = False #True # default is False; if True, compute dT decrements and return TESTING!!!!!!!!!!!; the only problem is that.... it doesn't compute covariance only stackedMap; we use this for anisotropic

mode = "kSZ"
if "tau_screening" in mode:
    filterType = "meanring" # og
    #filterType = "diskring" # TESTING!!!!!!!!!!!!
    Obs = 'tau'
elif "lensing" in mode:
    filterType = "meanring"
    Obs = 'tsz'
else:
    # filterType = "diskring"
    # filterType = 'DSigma'
    # filterType = stack_config.get("filterType", "DSigma")
    filterType = filterType # use the filter type from the config, which can be overridden by command line argument
    Obs = 'ksz' # og
    # Obs = 'ksz_anisotropic' # TESTING!!!! info hidden in vX, vY radian angle wrt theta
    #Obs = 'tsz_anisotropic' # TESTING!!!! info hidden in vX, vY radian angle wrt theta

# if want_sim:
#     pathMap = f'/pscratch/sd/b/boryanah/websky/abacus/map_tau_8192_ph201_fwhm1.6.fits'; extra_str += "_tau_Ill" #
#     #pathMap = f'/pscratch/sd/b/boryanah/websky/abacus/map_tau_8192_ph201_MTNG.fits'; extra_str += "_tau" #
#     #pathMap = f'/pscratch/sd/b/boryanah/websky/abacus/map_tau_8192_ph201_MTNG_spline.fits'; extra_str += "_tau_spline" # 0.45 - 0.55, halos at logM = 13.0
#     #pathMap = f'/pscratch/sd/b/boryanah/websky/abacus/lensed_map_tauvr_8192_ph201_MTNG.fits'; extra_str += "_tauvr"
#     #pathMap = f'/pscratch/sd/b/boryanah/websky/abacus/unlensed_map_tauvr_8192_ph201_MTNG.fits'; extra_str += "_unlensedtauvr"
#     #pathMap = f'/pscratch/sd/b/boryanah/websky/abacus/map_tauvr_8192_ph201_MTNG.fits'; extra_str += "_onlytauvr" #
#     #pathMap = f'/pscratch/sd/b/boryanah/websky/abacus/map_tauvr_8192_ph201_MTNG_spline.fits'; extra_str += "_onlytauvr_spline" # 0.45 - 0.55, halos at logM = 13.0
# else:
if stack_config.get('zeroV', False):
    print("Removing mean velocity from the catalog")
    extra_str += "_zeroV"

pathMap = "/pscratch/sd/b/boryanah/ACTxDESI/ACT/hilc_fullRes_TT_17000.fits" # 1.6 arcmin # OG
pathMask = '/pscratch/sd/b/boryanah/ACTxDESI/ACT/wide_mask_GAL070_apod_1.50_deg_wExtended_srcfree_Will.fits'
# TESTING: TODO: Maybe change back?
# pathMask = '/pscratch/sd/b/boryanah/ACTxDESI/ACT/wide_mask_GAL070_apod_1.50_deg_wExtended.fits'
# output_dir = f"/pscratch/sd/b/boryanah/ACTxDESI/output_test{extra_str}/"
# output_dir = f"/pscratch/sd/r/rhliu/projects/Weak_lensing/ksz_measurements/ACTxDESI/output_test{extra_str}/"
output_dir = output_dir_ + f"{extra_str}/"
pathHit = None

# catalogs = {"DESI_pz1": Catalog(u, massConversion, name="DESI_pz1", nameLong="DESI pz bin 1", out_dir=cat_dir, save=False, cat_fn=cat_fn)}
# catalogs = {"DESI_pz1": Catalog(u, massConversion, name="", nameLong="DESI pz bin 1", out_dir=cat_dir, save=True, cat_fn=cat_fn)}
# catalogs = {"DESI_pz1": Catalog(u, massConversion, name="DESI_pz1", nameLong="DESI pz bin 1", out_dir=cat_dir, save=False, cat_fn=cat_fn),
#             "DESI_pz2": Catalog(u, massConversion, name="DESI_pz2", nameLong="DESI pz bin 2", out_dir=cat_dir, save=False, cat_fn=cat_fn),
#             "DESI_pz3": Catalog(u, massConversion, name="DESI_pz3", nameLong="DESI pz bin 3", out_dir=cat_dir, save=False, cat_fn=cat_fn),
#             "DESI_pz4": Catalog(u, massConversion, name="DESI_pz4", nameLong="DESI pz bin 4", out_dir=cat_dir, save=False, cat_fn=cat_fn)}

catalogs = {}
for i, cat_name in enumerate(cat_names):
    catalogs[cat_name] = Catalog(u, massConversion, name=cat_name, nameLong=cat_name, out_dir=cat_dir, save=False, fig_dir=fig_dir, cat_fn=cat_fn)

# catalogs = {"DESIY3_LRG": Catalog(u, massConversion, name="DESIY3_LRG", nameLong="DESIY3 LRG all z bins", out_dir=cat_dir, save=False, fig_dir=fig_dir, cat_fn=cat_fn),
#             "DESIY3_z1": Catalog(u, massConversion, name="DESIY3_LRG_z1", nameLong="DESIY3 LRG z bin 1", out_dir=cat_dir, save=False, fig_dir=fig_dir, cat_fn=cat_fn),
#             "DESIY3_z2": Catalog(u, massConversion, name="DESIY3_LRG_z2", nameLong="DESIY3 LRG z bin 2", out_dir=cat_dir, save=False, fig_dir=fig_dir, cat_fn=cat_fn),
#             "DESIY3_z3": Catalog(u, massConversion, name="DESIY3_LRG_z3", nameLong="DESIY3 LRG z bin 3", out_dir=cat_dir, save=False, fig_dir=fig_dir, cat_fn=cat_fn),
#             "DESIY3_z4": Catalog(u, massConversion, name="DESIY3_LRG_z4", nameLong="DESIY3 LRG z bin 4", out_dir=cat_dir, save=False, fig_dir=fig_dir, cat_fn=cat_fn),
#             }
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

    #### CHANGE HERE IF NECESSARY: ZEROING VELOCITY COMPONENT:
    if stack_config.get('zeroV', False):
        print("Removing mean velocity from the catalog")
        catalog.vR -= np.mean(catalog.vR) # removing mean velocity # type: ignore
        # extra_str += "_zeroV"

    if stack_config.get('shuffle', False):
        print("Shuffling velocities in the catalog for null test")
        seed = stack_config.get('shuffle_seed', 42)
        rng = np.random.default_rng(seed=seed)
        rng.shuffle(catalog.vR) # shuffling velocity for null test # type: ignore
        
    # catalog.vR -= np.mean(catalog.vR) # removing mean velocity # type: ignore

    if Obs == 'tau':
        catalog.vR = catalog.vZ # hiding here info about T_large-scales

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
                    apod_pix=apod_pix)
    
    # Now we want to add some code to save the individual stacks
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