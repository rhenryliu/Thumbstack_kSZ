import sys
sys.path.append('/global/homes/b/boryanah/repos/ThumbStack')

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

##################################################################################

extra_str = "" # initialize
want_sim = False
if want_sim:
    #cat_dir = "/global/cfs/cdirs/desi/users/boryanah/kSZ_recon/for_fiona/"
    #cat_fn = "Extended_LRG_zerr0.0_AbacusSummit_huge_c000_ph201_masked.fits"
    #cat_fn = "Extended_LRG_zerr2.0_AbacusSummit_huge_c000_ph201_masked.fits"; extra_str += "_zerr2.0"
    #cat_fn = "Extended_LRG_zerr0.0_AbacusSummit_huge_c000_ph201_masked_tiny.fits"
    cat_dir = "/pscratch/sd/b/boryanah/websky/abacus/"
    cat_fn = "halos_LOGMhalo_Msunh13_masked.fits"
else:
    #cat_dir = "/global/cfs/cdirs/desi/users/boryanah/reconstruction_DESI/recon/"
    #cat_dir = '/pscratch/sd/b/boryanah/ACTxDESI/DESI/'
    cat_dir = '/pscratch/sd/b/boryanah/kSZ_pairwise/'
    #perc = 30
    perc = 10
    #cat_fn = "catalog_BGS_BRIGHT-20.2_R12.50_nmesh512_recsym_MG_masked.fits"
    #cat_fn = f"DESY6_ABSM{perc:d}.fits"
    #cat_fn = "ELG_LOPnotqso_masked.fits"
    #cat_fn = "ELG_LOPnotqso_cigale_masked.fits"
    #cat_fn = "LRG_cigale_masked.fits"
    cat_fn = "LRG_Y1_cigale_masked.fits"

# cosmological parameters
u = UnivMariana()
massConversion = MassConversionKravtsov14()
invPowerFunc = None
filterFuncRad = None
apod_pix = 20
save = True # if yes calculates delta T anew (saves mask, delta T, maybe vel is separate? but the numbers need to match)
#save = False # False means use saved stuff 
doBootstrap = True #False #True # do bootstrap or not? (not needed if e.g. shuffling velocities)
doMBins = False # we don't have mass bins yet
doVShuffle = False #True # default measurement is False
wantMF = False # default measurement is False
doOnlyFiltering = True #True # default is False; if True, compute dT decrements and return TESTING!!!!!!!!!!!; the only problem is that.... it doesn't compute covariance only stackedMap; we use this for anisotropic

mode = "kSZ"
if "tau_screening" in mode:
    filterType = "meanring" # og
    #filterType = "diskring" # TESTING!!!!!!!!!!!!
    Obs = 'tau'
elif "lensing" in mode:
    filterType = "meanring"
    Obs = 'tsz'
else:
    filterType = "diskring"
    #Obs = 'ksz' # og
    Obs = 'ksz_anisotropic' # TESTING!!!! info hidden in vX, vY radian angle wrt theta
    #Obs = 'tsz_anisotropic' # TESTING!!!! info hidden in vX, vY radian angle wrt theta

if want_sim:
    pathMap = f'/pscratch/sd/b/boryanah/websky/abacus/map_tau_8192_ph201_fwhm1.6.fits'; extra_str += "_tau_Ill" #
    #pathMap = f'/pscratch/sd/b/boryanah/websky/abacus/map_tau_8192_ph201_MTNG.fits'; extra_str += "_tau" #
    #pathMap = f'/pscratch/sd/b/boryanah/websky/abacus/map_tau_8192_ph201_MTNG_spline.fits'; extra_str += "_tau_spline" # 0.45 - 0.55, halos at logM = 13.0
    #pathMap = f'/pscratch/sd/b/boryanah/websky/abacus/lensed_map_tauvr_8192_ph201_MTNG.fits'; extra_str += "_tauvr"
    #pathMap = f'/pscratch/sd/b/boryanah/websky/abacus/unlensed_map_tauvr_8192_ph201_MTNG.fits'; extra_str += "_unlensedtauvr"
    #pathMap = f'/pscratch/sd/b/boryanah/websky/abacus/map_tauvr_8192_ph201_MTNG.fits'; extra_str += "_onlytauvr" # 
    #pathMap = f'/pscratch/sd/b/boryanah/websky/abacus/map_tauvr_8192_ph201_MTNG_spline.fits'; extra_str += "_onlytauvr_spline" # 0.45 - 0.55, halos at logM = 13.0
else:
    pathMap = "/pscratch/sd/b/boryanah/ACTxDESI/ACT/hilc_fullRes_TT_17000.fits" # 1.6 arcmin # OG
pathMask = '/pscratch/sd/b/boryanah/ACTxDESI/ACT/wide_mask_GAL070_apod_1.50_deg_wExtended_srcfree_Will.fits'
output_dir = f"/pscratch/sd/b/boryanah/ACTxDESI/output_test{extra_str}/"
pathHit = None

#catalogs = {"DESI_pz1": Catalog(u, massConversion, name="DESI_pz1", nameLong="DESI pz bin 1", out_dir=cat_dir, save=False, cat_fn=cat_fn)}
catalogs = {"DESI_pz1": Catalog(u, massConversion, name="", nameLong="DESI pz bin 1", out_dir=cat_dir, save=False, cat_fn=cat_fn)}

# Read CMB maps
nProc = 128

CMB_nu = 90.e9
cmap = cmbMap(pathMap,
              pathMask=pathMask,
              pathHit=pathHit,
              nu=CMB_nu, unitLatex=r'y',
              name='')

for key in catalogs.keys():
    catalog = catalogs[key]
    catalog.Mstellar = np.empty_like(catalog.RA)
    catalog.Mvir = np.empty_like(catalog.RA)
    catalog.integratedY = np.empty_like(catalog.RA)
    catalog.integratedKSZ = np.empty_like(catalog.RA)
    catalog.integratedTau = np.empty_like(catalog.RA)

    if Obs == 'tau':
        catalog.vR = catalog.vZ # hiding here info about T_large-scales

    ts = ThumbStack(u, catalog, 
                    cmap.map(),
                    cmap.mask(), 
                    cmap.hit(), 
                    '',
                    nameLong=None,
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

