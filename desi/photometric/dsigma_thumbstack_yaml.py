
import sys
sys.path.append('/global/homes/r/rhliu/projects/repos/ThumbStack')

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
# param_Dict = json.loads(sys.argv[1])

##################################################################################
t0 = time.time()

extra_str = "" # initialize
want_sim = False

parser = argparse.ArgumentParser(description='Process config.')
parser.add_argument('-p', '--path2config', type=str, default='./configs/diskring_NGC.yaml', help='Path to the configuration file.')
args = vars(parser.parse_args())
print(f"Arguments: {args}")
path2config = args['path2config']
with open(path2config, 'r') as f:
    config = yaml.safe_load(f)
    
cat_dir = config.get('cat_dir', '/pscratch/sd/r/rhliu/projects/Weak_lensing/desi/spec_Y3/')
fig_dir = config.get('fig_dir', '/pscratch/sd/r/rhliu/projects/Weak_lensing/figs/')
cat_fn_template = config.get('cat_fn', 'catalog_{field}_{filter_cut}.txt')
field = config.get('field', 'NGC') # NGC, SGC, or full
filter_cut = config.get('filter_cut', 'unfiltered') # unfiltered or no_src_with_cluster_mask
cat_fn = cat_fn_template.format(field=field, filter_cut=filter_cut)
extra_str += config.get('extra_str', f"/{field}_{filter_cut}").format(field=field, filter_cut=filter_cut)

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
save = config.get('save', False) # if yes calculates delta T anew (saves mask, delta T, maybe vel is separate? but the numbers need to match)
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
    filterType = config.get("filterType", "DSigma")
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
if config.get('zeroV', False):
    print("Removing mean velocity from the catalog")
    extra_str += "_zeroV"

pathMap = "/pscratch/sd/b/boryanah/ACTxDESI/ACT/hilc_fullRes_TT_17000.fits" # 1.6 arcmin # OG
pathMask = '/pscratch/sd/b/boryanah/ACTxDESI/ACT/wide_mask_GAL070_apod_1.50_deg_wExtended_srcfree_Will.fits'
# output_dir = f"/pscratch/sd/b/boryanah/ACTxDESI/output_test{extra_str}/"
# output_dir = f"/pscratch/sd/r/rhliu/projects/Weak_lensing/ksz_measurements/ACTxDESI/output_test{extra_str}/"
output_dir = f"/pscratch/sd/r/rhliu/projects/Weak_lensing/ksz_measurements/ACTxDESI/spec_Y3{extra_str}/"
pathHit = None

# catalogs = {"DESI_pz1": Catalog(u, massConversion, name="DESI_pz1", nameLong="DESI pz bin 1", out_dir=cat_dir, save=False, cat_fn=cat_fn)}
# catalogs = {"DESI_pz1": Catalog(u, massConversion, name="", nameLong="DESI pz bin 1", out_dir=cat_dir, save=True, cat_fn=cat_fn)}
# catalogs = {"DESI_pz1": Catalog(u, massConversion, name="DESI_pz1", nameLong="DESI pz bin 1", out_dir=cat_dir, save=False, cat_fn=cat_fn),
#             "DESI_pz2": Catalog(u, massConversion, name="DESI_pz2", nameLong="DESI pz bin 2", out_dir=cat_dir, save=False, cat_fn=cat_fn),
#             "DESI_pz3": Catalog(u, massConversion, name="DESI_pz3", nameLong="DESI pz bin 3", out_dir=cat_dir, save=False, cat_fn=cat_fn),
#             "DESI_pz4": Catalog(u, massConversion, name="DESI_pz4", nameLong="DESI pz bin 4", out_dir=cat_dir, save=False, cat_fn=cat_fn)}

catalogs = {"DESIY3_LRG": Catalog(u, massConversion, name="DESIY3_LRG", nameLong="DESIY3 LRG all z bins", out_dir=cat_dir, save=False, fig_dir=fig_dir, cat_fn=cat_fn),
            # "DESIY3_z1": Catalog(u, massConversion, name="DESIY3_LRG_z1", nameLong="DESIY3 LRG z bin 1", out_dir=cat_dir, save=False, fig_dir=fig_dir, cat_fn=cat_fn),
            # "DESIY3_z2": Catalog(u, massConversion, name="DESIY3_LRG_z2", nameLong="DESIY3 LRG z bin 2", out_dir=cat_dir, save=False, fig_dir=fig_dir, cat_fn=cat_fn),
            # "DESIY3_z3": Catalog(u, massConversion, name="DESIY3_LRG_z3", nameLong="DESIY3 LRG z bin 3", out_dir=cat_dir, save=False, fig_dir=fig_dir, cat_fn=cat_fn),
            # "DESIY3_z4": Catalog(u, massConversion, name="DESIY3_LRG_z4", nameLong="DESIY3 LRG z bin 4", out_dir=cat_dir, save=False, fig_dir=fig_dir, cat_fn=cat_fn),
            }
# Read CMB maps
nProc = 128

CMB_nu = 90.e9
cmap = cmbMap(pathMap,
              pathMask=pathMask,
              pathHit=pathHit,
              nu=CMB_nu, unitLatex=r'y',
              name='')

for i, key in enumerate(catalogs.keys()):
    catalog = catalogs[key]
    catalog.Mstellar = np.empty_like(catalog.RA)
    catalog.Mvir = np.empty_like(catalog.RA)
    catalog.integratedY = np.empty_like(catalog.RA)
    catalog.integratedKSZ = np.empty_like(catalog.RA)
    catalog.integratedTau = np.empty_like(catalog.RA)

    #### CHANGE HERE IF NECESSSARY: ZEROING VELOCITY COMPONENT:
    if config.get('zeroV', False):
        # print("Removing mean velocity from the catalog")
        catalog.vR -= np.mean(catalog.vR) # removing mean velocity # type: ignore
        # extra_str += "_zeroV"
        
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
                    output_dir=output_dir,
                    Obs=Obs,
                    wantMF=wantMF,
                    invPowerFunc=invPowerFunc,
                    filterFuncRad=filterFuncRad,
                    apod_pix=apod_pix)
    
    # Now we want to add some code to save the individual stacks
    if config.get('save_invidiual_stacks', False):
        cat_fn_dataframe_template = config.get('cat_fn_dataframe', '{field}_zbin{bin}_{filter_cut}.csv')
        # df_catalog = pd.read_csv(cat_dir + cat_fn, delim_whitespace=True)
        data_names = ['TARGETID', 'Z','RA', 'DEC', 'VEL_LOS_RENORM']
        data_df = pd.read_csv(cat_dir + cat_fn_dataframe_template.format(field=field, bin=i+1, filter_cut=filter_cut))
        data_df = data_df[data_names]
        data_df.rename(columns={'VEL_LOS_RENORM': 'vR'}, inplace=True) # type: ignore

        allProfiles = computeProfiles(ts, filterType, est='ksz_uniformweight')
        print("Computed all profiles for individual stacks. Allprofiles shape:", allProfiles.shape)
        mask = ts.catalogMask(overlap=True, psMask=True, filterType=filterType, mVir=None)
        
        # data_df['Mask'] = mask.astype(int)
        data_df = data_df[mask]
        NObj = data_df.shape[0]

        for j, R in enumerate(ts.RApArcmin):
            print(str(R), )
            stacks = np.zeros(NObj)
            # stacks[mask] = allProfiles[:, j]
            stacks = allProfiles[:, j]

            data_df[str(R)] = stacks

        save_name = config.get('save_dir', '').format(field=field, filter_cut=filter_cut) + config.get('save_name', '{filterType}_profiles_{field}_{filter_cut}_z{bin}.csv').format(filterType=filterType, field=field, bin=i+1, filter_cut=filter_cut)
        print(f"Saving individual stacks to {save_name}")
        
        data_df.to_csv(save_name, index=False)
        
t1 = time.time()
print(f"Total time: {t1 - t0} seconds")