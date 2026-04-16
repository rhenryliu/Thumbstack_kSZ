import pandas as pd
import numpy as np
from astropy.table import Table
from pixell import enmap
from cosmoprimo.fiducial import DESI # type: ignore
import pyclass # type: ignore
import yaml
import argparse
import healpy
import os
import sys

sys.path.append('./')
from remove_close_pairs import enforce_min_separation, physical_to_min_sep, mean_redshift_from_table
import astropy.units as u
from astropy.table import Table as AstropyTable

# ---------------------------------------------------------------------------
# Configuration
# ---------------------------------------------------------------------------

parser = argparse.ArgumentParser(description='Prepare DES/non-DES split SGC catalogues from DESI data.')
parser.add_argument('-p', '--path2config', type=str,
                    default='./configs/DES/prepare_cat_LRG_zeroV_DES.yaml',
                    help='Path to the configuration file.')
args = vars(parser.parse_args())
print(f"Arguments: {args}")
path2config = args['path2config']
with open(path2config, 'r') as f:
    config = yaml.safe_load(f)

main_directory = config['source'].get('main_directory', '/global/cfs/cdirs/desi/survey/catalogs/DA2/')
post_rec_directory = config['source'].get('post_rec_dir', 'analysis/loa-v1/LSScats/v1.1/BAO/unblinded/desipipe/2pt/recon_sm15_IFFT_recsym/')
pre_rec_directory = config['source'].get('pre_rec_dir', 'LSS/loa-v1/LSScats/v1.1/nonKP/')
NGC_fn = config['source'].get('NGC_fn', 'LRG_NGC_clustering.dat.fits')
SGC_fn = config['source'].get('SGC_fn', 'LRG_SGC_clustering.dat.fits')

cat_type = config['processing'].get('cat_type', 'LRG')
filter_type = config['processing'].get('filter_type', 'no_src_with_cluster_mask')
zeroV = config['processing'].get('zeroV', False)
zbins_config = config['processing'].get('zbins', None)
if zbins_config is not None:
    zbins = [tuple(z) for z in config["processing"]["zbins"]]
else:
    zbins = None

# Fields to produce. 'SGC_DES' and 'SGC_nonDES' are always the core outputs.
# Add 'NGC' and/or 'full' to also produce those (Interpretation A: full = NGC + SGC).
# 'full' requires 'NGC' to be present as well.
fields = config['processing'].get('fields', ['SGC_DES', 'SGC_nonDES'])

output_dir = config['save'].get('save_dir', '/pscratch/sd/r/rhliu/projects/Weak_lensing/desi/spec_Y3/LRG_catalogues_zeroV/')
save_format = config['save'].get('save_format', 'csv')

apply_act_filter = filter_type in ['no_src_with_cluster_mask', 'nopairs']
apply_pair_pruning = filter_type == 'nopairs'
if apply_pair_pruning:
    physical_sep = config['processing'].get('physical_sep', 1.0)

need_NGC = ('NGC' in fields) or ('full' in fields) or ('full_nonDES' in fields)

print(f"Fields to produce: {fields}")
if zeroV:
    print("zeroV=True: mean subtraction applied per field and per z-bin independently.")
else:
    print("Using original velocities without zeroV correction.")

# ---------------------------------------------------------------------------
# Load DESI FITS catalogues
# ---------------------------------------------------------------------------

# SGC: always loaded (needed for DES mask split)
print("\nLoading SGC catalogues...")
dat_pre_rec_SGC = Table.read(main_directory + pre_rec_directory + SGC_fn, format='fits')
pre_rec_SGC = dat_pre_rec_SGC.to_pandas()
dat_post_rec_SGC = Table.read(main_directory + post_rec_directory + SGC_fn, format='fits')
post_rec_SGC = dat_post_rec_SGC.to_pandas()
pre_rec_SGC["GC"] = "SGC"

# NGC: only loaded when 'NGC' or 'full' is requested
if need_NGC:
    print("Loading NGC catalogues...")
    dat_pre_rec_NGC = Table.read(main_directory + pre_rec_directory + NGC_fn, format='fits')
    pre_rec_NGC = dat_pre_rec_NGC.to_pandas()
    dat_post_rec_NGC = Table.read(main_directory + post_rec_directory + NGC_fn, format='fits')
    post_rec_NGC = dat_post_rec_NGC.to_pandas()
    pre_rec_NGC["GC"] = "NGC"

# ---------------------------------------------------------------------------
# Velocity computation
# ---------------------------------------------------------------------------

cosmo = DESI()
z_eff = 0.780
f_eff = cosmo.growth_rate(z_eff)

print("\nComputing SGC velocities...")
chi_pre_SGC = cosmo.comoving_radial_distance(pre_rec_SGC['Z'])   # Mpc/h
chi_post_SGC = cosmo.comoving_radial_distance(post_rec_SGC['Z']) # Mpc/h
disp_SGC_LOS = -(chi_post_SGC - chi_pre_SGC) / cosmo.h          # Mpc
disp_SGC_LOS /= (1 + f_eff)                                      # Mpc, RSD correction
pre_rec_SGC["DISP_LOS"] = disp_SGC_LOS

f_SGC = cosmo.growth_rate(post_rec_SGC['Z'])
H_z_SGC = cosmo.hubble_function(post_rec_SGC['Z'])  # km/(s Mpc)
a_SGC = 1 / (1 + post_rec_SGC['Z'])
vel_SGC_LOS = a_SGC * H_z_SGC * f_SGC * disp_SGC_LOS  # km/s
pre_rec_SGC["VEL_LOS"] = vel_SGC_LOS
pre_rec_SGC["VEL_LOS_RENORM"] = vel_SGC_LOS * cosmo.growth_factor(pre_rec_SGC["Z"]) / cosmo.growth_factor(z_eff)

if need_NGC:
    print("Computing NGC velocities...")
    chi_pre_NGC = cosmo.comoving_radial_distance(pre_rec_NGC['Z'])
    chi_post_NGC = cosmo.comoving_radial_distance(post_rec_NGC['Z'])
    disp_NGC_LOS = -(chi_post_NGC - chi_pre_NGC) / cosmo.h
    disp_NGC_LOS /= (1 + f_eff)
    pre_rec_NGC["DISP_LOS"] = disp_NGC_LOS

    f_NGC = cosmo.growth_rate(post_rec_NGC['Z'])
    H_z_NGC = cosmo.hubble_function(post_rec_NGC['Z'])
    a_NGC = 1 / (1 + post_rec_NGC['Z'])
    vel_NGC_LOS = a_NGC * H_z_NGC * f_NGC * disp_NGC_LOS
    pre_rec_NGC["VEL_LOS"] = vel_NGC_LOS
    pre_rec_NGC["VEL_LOS_RENORM"] = vel_NGC_LOS * cosmo.growth_factor(pre_rec_NGC["Z"]) / cosmo.growth_factor(z_eff)

# Sort by redshift
pre_rec_SGC_sort = pd.DataFrame(pre_rec_SGC).sort_values("Z")
if need_NGC:
    pre_rec_NGC_sort = pd.DataFrame(pre_rec_NGC).sort_values("Z")

# ---------------------------------------------------------------------------
# ACT mask filtering
# ---------------------------------------------------------------------------

def sky2map(ra, dec, cmbMap):
    """Map value at (ra, dec) in degrees via nearest-neighbor interpolation."""
    sourcecoord = np.array([dec, ra]) * (np.pi / 180)
    return cmbMap.at(sourcecoord, order=0)

def apply_act_overlap_filter(catalog_df, cmbMask, thresh=0.95, name="catalog"):
    ra = np.array(catalog_df["RA"])
    dec = np.array(catalog_df["DEC"])
    nObj = len(ra)
    print(f"Applying ACT overlap filter to {name}: {nObj:,} objects")
    hit = sky2map(ra, dec, cmbMask)
    overlapFlag = np.array(hit > thresh) * 1
    filtered = catalog_df[overlapFlag == 1]
    print(f"After ACT filtering - {name}: {len(filtered):,} objects ({len(filtered)/nObj*100:.1f}%)")
    return filtered

if apply_act_filter:
    print("\n=== Applying ACT overlap filtering ===")
    masks_directory = '/pscratch/sd/j/jia_qu/ACTxDESIY3/'
    cmbMask = enmap.read_fits(f"{masks_directory}/wide_mask_GAL070_apod_1.50_deg_wExtended_no_src_with_cluster.fits")
    pre_rec_SGC_sort_ACT = apply_act_overlap_filter(pre_rec_SGC_sort, cmbMask, name="SGC")
    if need_NGC:
        pre_rec_NGC_sort_ACT = apply_act_overlap_filter(pre_rec_NGC_sort, cmbMask, name="NGC")
else:
    print("\n=== Skipping ACT overlap filtering ===")
    pre_rec_SGC_sort_ACT = pre_rec_SGC_sort
    if need_NGC:
        pre_rec_NGC_sort_ACT = pre_rec_NGC_sort

# ---------------------------------------------------------------------------
# DES footprint mask — split SGC into DES and non-DES subregions
# ---------------------------------------------------------------------------

print("\n=== Applying DES footprint mask to SGC ===")
# Mask is binary (0 = outside DES footprint, 1 = inside), NSIDE=256 HEALPix
DESMask = healpy.fitsfunc.read_map('./masks/mask_desy3_footprint.fits')
nside_des = healpy.npix2nside(len(DESMask))

ra_sgc = np.array(pre_rec_SGC_sort_ACT["RA"])
dec_sgc = np.array(pre_rec_SGC_sort_ACT["DEC"])
pix_sgc = healpy.ang2pix(nside_des, ra_sgc, dec_sgc, lonlat=True)
in_des = DESMask[pix_sgc] > 0

pre_rec_SGC_DES = pre_rec_SGC_sort_ACT[in_des].copy()
pre_rec_SGC_nonDES = pre_rec_SGC_sort_ACT[~in_des].copy()
n_sgc = len(pre_rec_SGC_sort_ACT)
print(f"SGC in DES footprint    (SGC_DES):    {len(pre_rec_SGC_DES):,} objects ({len(pre_rec_SGC_DES)/n_sgc*100:.1f}%)")
print(f"SGC outside DES footprint (SGC_nonDES): {len(pre_rec_SGC_nonDES):,} objects ({len(pre_rec_SGC_nonDES)/n_sgc*100:.1f}%)")

# ---------------------------------------------------------------------------
# Internal primary caps — independently zeroed
# SGC_DES and SGC_nonDES are always computed (complementary pair from DES mask split).
# NGC is computed whenever any assembled field needs it ('NGC', 'full', 'full_nonDES').
# Each cap is independently zeroed; derived fields (full, full_nonDES) are assembled
# from the already-zeroed caps and are not re-zeroed.
# ---------------------------------------------------------------------------

_internal_caps = {'SGC_DES': pre_rec_SGC_DES, 'SGC_nonDES': pre_rec_SGC_nonDES}
if need_NGC:
    _internal_caps['NGC'] = pre_rec_NGC_sort_ACT
    if 'NGC' not in fields:
        print("(NGC loaded internally for assembled fields; not saved as a standalone output)")

# Non-binned zeroV: subtract mean VEL_LOS_RENORM independently for each cap.
if zeroV:
    print("\n=== Applying non-binned zeroV corrections ===")
    for cap in list(_internal_caps.keys()):
        cat = _internal_caps[cap]
        if len(cat) > 0:
            _internal_caps[cap] = cat.copy()
            mean_v = np.mean(_internal_caps[cap]["VEL_LOS_RENORM"])
            _internal_caps[cap]["VEL_LOS_RENORM"] -= mean_v
            print(f"  {cap}: subtracted mean={mean_v:.2f} km/s from {len(_internal_caps[cap]):,} objects")

# Assemble output catalog dict from zeroed caps.
# Primary caps that are output fields pass through directly.
# Derived fields (full, full_nonDES) are concatenations of zeroed primaries.
all_catalogs = {}
for cap in ['SGC_DES', 'SGC_nonDES', 'NGC']:
    if cap in fields and cap in _internal_caps:
        all_catalogs[cap] = _internal_caps[cap]

if 'full' in fields and need_NGC:
    components = [_internal_caps[c] for c in ['NGC', 'SGC_DES', 'SGC_nonDES'] if c in _internal_caps]
    all_catalogs['full'] = pd.concat(components).sort_values("Z")
    component_names = [c for c in ['NGC', 'SGC_DES', 'SGC_nonDES'] if c in _internal_caps]
    print(f"  full: assembled from {', '.join(component_names)} ({len(all_catalogs['full']):,} objects)")

if 'full_nonDES' in fields and need_NGC:
    components = [_internal_caps[c] for c in ['NGC', 'SGC_nonDES'] if c in _internal_caps]
    all_catalogs['full_nonDES'] = pd.concat(components).sort_values("Z")
    component_names = [c for c in ['NGC', 'SGC_nonDES'] if c in _internal_caps]
    print(f"  full_nonDES: assembled from {', '.join(component_names)} ({len(all_catalogs['full_nonDES']):,} objects)")

# ---------------------------------------------------------------------------
# Redshift bins
# ---------------------------------------------------------------------------

def create_redshift_bins(catalog_df, z_bins, name):
    binned = {}
    for i, (z_min, z_max) in enumerate(z_bins, 1):
        bin_cat = catalog_df[(catalog_df['Z'] > z_min) & (catalog_df['Z'] <= z_max)]
        binned[f"{name}_zbin{i}"] = bin_cat
        print(f"  Z-bin {i} ({z_min:.2f}-{z_max:.2f}): {len(bin_cat):,} objects")
    return binned

print("\n=== Creating redshift bins ===")
_internal_zbins = {}
for cap, cat in _internal_caps.items():
    print(f"\n{cap}:")
    zbins_cap = create_redshift_bins(cat, zbins, cap)
    if zeroV:
        # Per-bin zeroV: subtract mean within each bin independently
        for key in zbins_cap:
            bin_cat = zbins_cap[key]
            if len(bin_cat) > 0:
                zbins_cap[key] = bin_cat.copy()
                mean_v = np.mean(zbins_cap[key]["VEL_LOS_RENORM"])
                zbins_cap[key]["VEL_LOS_RENORM"] -= mean_v
    _internal_zbins[cap] = zbins_cap

# Assemble output z-bins from zeroed internal cap z-bins.
all_zbins = {}
for cap in ['SGC_DES', 'SGC_nonDES', 'NGC']:
    if cap in fields and cap in _internal_zbins:
        all_zbins[cap] = _internal_zbins[cap]

n_zbins = len(zbins) # type: ignore

if 'full' in fields and need_NGC:
    print(f"\nfull (assembled from component z-bins):")
    full_zbins = {}
    for i in range(1, n_zbins + 1):
        components = [_internal_zbins[c][f"{c}_zbin{i}"] for c in ['NGC', 'SGC_DES', 'SGC_nonDES'] if c in _internal_zbins]
        full_zbin = pd.concat(components).sort_values("Z")
        full_zbins[f"full_zbin{i}"] = full_zbin
        z_min, z_max = zbins[i - 1] # type: ignore
        print(f"  Z-bin {i} ({z_min:.2f}-{z_max:.2f}): {len(full_zbin):,} objects")
    all_zbins['full'] = full_zbins

if 'full_nonDES' in fields and need_NGC:
    print(f"\nfull_nonDES (assembled from component z-bins):")
    full_nonDES_zbins = {}
    for i in range(1, n_zbins + 1):
        components = [_internal_zbins[c][f"{c}_zbin{i}"] for c in ['NGC', 'SGC_nonDES'] if c in _internal_zbins]
        full_nonDES_zbin = pd.concat(components).sort_values("Z")
        full_nonDES_zbins[f"full_nonDES_zbin{i}"] = full_nonDES_zbin
        z_min, z_max = zbins[i - 1] # type: ignore
        print(f"  Z-bin {i} ({z_min:.2f}-{z_max:.2f}): {len(full_nonDES_zbin):,} objects")
    all_zbins['full_nonDES'] = full_nonDES_zbins

# ---------------------------------------------------------------------------
# Save catalogs
# ---------------------------------------------------------------------------

def save_catalogs():
    os.makedirs(output_dir, exist_ok=True)

    if apply_pair_pruning:
        print(f"\n=== Pruning and saving to {output_dir} ===")
        print(f"Physical separation threshold: {physical_sep} Mpc")
        # Conservative z_mean: lowest z-bin of the first internal cap
        first_cap = next(iter(_internal_caps))
        first_zbin = _internal_zbins[first_cap][f"{first_cap}_zbin1"]
        z_mean_conservative = mean_redshift_from_table(AstropyTable.from_pandas(first_zbin), z_col="Z")
        print(f"Using conservative z_mean = {z_mean_conservative:.4f} for full-z catalogs")
    else:
        print(f"\n=== Saving catalogs to {output_dir} (no pruning) ===")

    for fname in all_catalogs.keys():
        catalogs_to_save = {f"{fname}_catalog_Y3": all_catalogs[fname]}
        catalogs_to_save.update(all_zbins[fname])

        for name, catalog in catalogs_to_save.items():
            if len(catalog) == 0:
                print(f"Skipping {name}: empty catalog")
                continue

            if apply_pair_pruning:
                cat_table = AstropyTable.from_pandas(catalog)
                z_mean = mean_redshift_from_table(cat_table, z_col="Z")
                if "zbin" not in name:
                    z_mean = z_mean_conservative
                min_sep = physical_to_min_sep(physical_sep * u.Mpc, z_mean, comoving=False, out_unit=u.arcmin)
                print(f"\n{name}: z_mean={z_mean:.4f}, min_sep={min_sep:.3f}")
                print(f"  Before pruning: {len(cat_table):,} objects")
                pruned_table = enforce_min_separation(cat_table, ra_col="RA", dec_col="DEC", min_sep=min_sep)
                print(f"  After pruning: {len(pruned_table):,} objects ({len(pruned_table)/len(cat_table)*100:.1f}%)")
                final_catalog = pruned_table.to_pandas() # type: ignore
            else:
                final_catalog = catalog

            output_file = f"{output_dir}/{name}_{filter_type}.{save_format}"
            if save_format == "csv":
                final_catalog.to_csv(output_file, index=False)
            elif save_format == "txt":
                np.savetxt(output_file, np.array(final_catalog))
            print(f"Saved {name}: {len(final_catalog):,} objects -> {output_file}")

save_catalogs()

print(f"\n=== Summary ===")
for fname, cat in all_catalogs.items():
    total_bins = sum(len(c) for c in all_zbins[fname].values())
    print(f"{fname}: {len(cat):,} objects | z-bins total: {total_bins:,}")
print('done!!')
