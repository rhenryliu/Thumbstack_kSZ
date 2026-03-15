import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from astropy.table import Table
from joblib import Parallel, delayed
from tqdm import tqdm
import time
from pixell import enmap
import copy
import matplotlib
from scipy import special, optimize, integrate, stats
# import classy module
from classy import Class # type: ignore
from scipy import integrate
from cosmoprimo.fiducial import DESI # tienes que tener el environment de cosmodesi # type: ignore
import pandas as pd
import pyclass # type: ignore
import yaml
import argparse

# Import pruning functions from remove_close_pairs
import sys
sys.path.append('./')
from remove_close_pairs import enforce_min_separation, physical_to_min_sep, mean_redshift_from_table
import astropy.units as u
from astropy.table import Table as AstropyTable

# Load configuration
parser = argparse.ArgumentParser(description='Process config.')
parser.add_argument('-p', '--path2config', type=str, default='./configs/prepare_cat_BGS_cigale.yaml', help='Path to the configuration file.')
args = vars(parser.parse_args())
print(f"Arguments: {args}")
path2config = args['path2config']
with open(path2config, 'r') as f:
    config = yaml.safe_load(f)

# Extract configuration values
cat_fn = config['source'].get('cat_fn', '')

cat_type = config['processing'].get('cat_type', 'LRG')
filter_type = config['processing'].get('filter_type', 'nopairs')
zbins_config = config['processing'].get('zbins', None)
if zbins_config is not None:
    zbins = [tuple(z) for z in config["processing"]["zbins"]]
else:
    zbins = None

output_dir = config['save'].get('save_dir', '/pscratch/sd/r/rhliu/projects/Weak_lensing/desi/spec_Y3/')
save_format = config['save'].get('save_format', 'csv')

# Load single catalogue
dat = Table.read(cat_fn, format='fits')
pre_rec = dat.to_pandas()

# vR already contains the reconstructed LOS velocity; expose it as VEL_LOS_RENORM
pre_rec["VEL_LOS_RENORM"] = pre_rec["vR"]

# We sort them by redshift
pre_rec_sort = pd.DataFrame(pre_rec).sort_values("Z")


# Select the ones overlapping with ACT

def sky2map(ra, dec, cmbMap):
    '''Gives the map value at coordinates (ra, dec).
    ra, dec in degrees.
    Uses nearest neighbor, no interpolation.
    Will return 0 if the coordinates requested are outside the map
    '''
    # interpolate the map to the given sky coordinates
    sourcecoord = np.array([dec, ra]) * (np.pi / 180)   # convert from degrees to radians
    # use nearest neighbor interpolation
    return cmbMap.at(sourcecoord, order=0)

def apply_act_overlap_filter(catalog_df, cmbMask, thresh=0.95, name="catalog"):
    """
    Apply ACT overlap filtering to a catalog DataFrame
    
    Parameters:
    -----------
    catalog_df : pandas.DataFrame
        Input catalog with RA, DEC columns
    cmbMask : enmap
        CMB mask map
    thresh : float
        Threshold for overlap (default 0.95)
    name : str
        Name for logging purposes
        
    Returns:
    --------
    pandas.DataFrame : Filtered catalog
    """
    ra = np.array(catalog_df["RA"])
    dec = np.array(catalog_df["DEC"])
    nObj = len(ra)
    
    print(f"Applying ACT overlap filter to {name}: {nObj} objects")
    
    # Vectorized approach for better performance
    hit = sky2map(ra, dec, cmbMask)
    overlapFlag = np.array(hit>thresh)*1
    
    filtered_catalog = catalog_df[overlapFlag==1]
    print(f"After ACT filtering - {name}: {len(filtered_catalog)} objects ({len(filtered_catalog)/nObj*100:.1f}%)")
    
    return filtered_catalog

def create_redshift_bins(catalog_df, z_bins=None, name="catalog"):
    """
    Create redshift bins from a catalog
    
    Parameters:
    -----------
    catalog_df : pandas.DataFrame
        Input catalog with Z column
    z_bins : list of tuples
        List of (z_min, z_max) tuples. If None, uses default bins
    name : str
        Name for logging purposes
        
    Returns:
    --------
    dict : Dictionary with bin names as keys and DataFrames as values
    """
    if z_bins is None:
        z_bins = [
            (0.4, 0.6),
            (0.6, 0.8), 
            (0.8, 0.95),
            (0.95, 1.1)
        ]
    
    binned_catalogs = {}
    
    for i, (z_min, z_max) in enumerate(z_bins, 1):
        bin_catalog = catalog_df[(catalog_df['Z'] > z_min) & (catalog_df['Z'] <= z_max)]
        bin_name = f"{name}_zbin{i}"
        binned_catalogs[bin_name] = bin_catalog
        print(f"Z-bin {i} ({z_min}-{z_max}): {len(bin_catalog)} objects")
    
    return binned_catalogs

# Load CMB maps
print("Loading CMB maps...")
masks_directory = '/pscratch/sd/j/jia_qu/ACTxDESIY3/'
cmbMap = enmap.read_fits(f"{masks_directory}/hilc_fullRes_TT_17000.fits")
#cmbMask = enmap.read_fits("/home/jiaqu/Thumbstack_DESI/wide_mask_GAL070_apod_1.50_deg_wExtended_srcfree_Will.fits")
#cmbMask = enmap.read_fits("/project/rrg-rbond-ac/msyriac/ilc_dr6v3/20230606/wide_mask_GAL070_apod_1.50_deg_wExtended.fits")
cmbMask = enmap.read_fits(f"{masks_directory}/wide_mask_GAL070_apod_1.50_deg_wExtended_no_src_with_cluster.fits")

# Apply ACT overlap filtering based on filter_type
print(f"\n=== Processing with filter_type: {filter_type} ===")
apply_act_filter = filter_type in ['no_src_with_cluster_mask', 'nopairs']
apply_pair_pruning = filter_type == 'nopairs'
if apply_pair_pruning:
    physical_sep = config['processing'].get('physical_sep', 1.0) # Mpc

if apply_act_filter:
    print("\n=== Applying ACT Overlap Filtering ===")
    pre_rec_sort_ACT = apply_act_overlap_filter(pre_rec_sort, cmbMask, name="full catalog")
else:
    print("\n=== Skipping ACT Overlap Filtering ===")
    pre_rec_sort_ACT = pre_rec_sort

print(pre_rec_sort_ACT.columns)
print("Total objects in full catalog:")
print(len(pre_rec_sort_ACT))

# Create redshift bins for all filtered catalogs
print("\n=== Creating Redshift Bins ===")
print("Full catalog bins:")
full_zbins = create_redshift_bins(pre_rec_sort_ACT, z_bins=zbins, name="full")

# For backward compatibility, keep the original variable names
df = pre_rec_sort_ACT

def save_catalogs():
    """
    Save all processed catalogs to files
    
    """
    import os
    os.makedirs(output_dir, exist_ok=True)
    
    # Save main catalogs
    catalogs_to_save = {
        'full_catalog_Y3': pre_rec_sort_ACT,
    }
    
    # Add redshift bins
    catalogs_to_save.update(full_zbins)
    
    if apply_pair_pruning:
        print(f"\n=== Pruning and Saving Catalogs to {output_dir} ===")
        print(f"Physical separation threshold: {physical_sep} Mpc")
        
        # Compute mean redshift of lowest z-bin for conservative full catalog pruning
        lowest_zbin = full_zbins['full_zbin1']
        z_mean_conservative = mean_redshift_from_table(AstropyTable.from_pandas(lowest_zbin), z_col="Z")
        print(f"Using conservative z_mean = {z_mean_conservative:.4f} (from lowest z-bin) for full catalogs")
    else:
        print(f"\n=== Saving Catalogs to {output_dir} (no pruning) ===")
    
    for name, catalog in catalogs_to_save.items():
        if len(catalog) > 0:
            if apply_pair_pruning:
                # Convert to astropy Table for pruning
                cat_table = AstropyTable.from_pandas(catalog)
                
                # Compute mean redshift for this catalog
                z_mean = mean_redshift_from_table(cat_table, z_col="Z")
                
                # Use conservative z_mean for full catalogs (not z-bins)
                if "zbin" not in name:
                    z_mean = z_mean_conservative
                
                # Convert physical separation to angular separation at z_mean
                min_sep = physical_to_min_sep(physical_sep * u.Mpc, z_mean, comoving=False, out_unit=u.arcmin)
                
                print(f"\n{name}: z_mean={z_mean:.4f}, min_sep={min_sep:.3f}")
                print(f"  Before pruning: {len(cat_table)} objects")
                
                # Prune catalog
                pruned_table = enforce_min_separation(
                    cat_table, 
                    ra_col="RA", 
                    dec_col="DEC", 
                    min_sep=min_sep
                )
                
                print(f"  After pruning: {len(pruned_table)} objects ({len(pruned_table)/len(cat_table)*100:.1f}%)")
                
                # Convert back to pandas
                final_catalog = pruned_table.to_pandas() # type: ignore
            else:
                final_catalog = catalog
            
            # Construct output filename based on filter_type
            output_file = f"{output_dir}/{name}_{filter_type}.{save_format}"
            
            if save_format == "csv":
                final_catalog.to_csv(output_file, index=False)
                print(f"Saved {name}: {len(final_catalog)} objects -> {output_file}")
            elif save_format == "txt":
                np.savetxt(output_file, np.array(final_catalog))
                print(f"Saved {name}: {len(final_catalog)} objects -> {output_file}")
        else:
            print(f"Skipping {name}: empty catalog")

# Save catalogs
save_catalogs()

print(f"\n=== Summary ===")
print(f"Full catalog: {len(pre_rec_sort_ACT)} objects")
print(f"Total in z-bins: {sum(len(cat) for cat in full_zbins.values())} objects")
print('done!!')
