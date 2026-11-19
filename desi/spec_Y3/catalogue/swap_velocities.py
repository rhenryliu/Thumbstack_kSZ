"""Swap the velocities of the production LRG (or BGS) catalogues for another velocity estimate.

The production per-bin catalogues ({field}_zbin{i}_{filter_type}.csv, written by prepare_catalogue_yaml.py)
are kept object for object and in the same row order, and only VEL_LOS and VEL_LOS_RENORM are replaced, so
the stacks can reuse the production per-object filtering (stack.filtering_dir in make_stacks_yaml.py).
convert_to_TS_catalog.py then writes the ThumbStack catalogues from the same config.

velocity_column:
    'official': the production VEL_LOS (BAO RecSym reconstruction), the known-answer control: with
        renorm and zeroV it must reproduce the production VEL_LOS_RENORM.
    any column of velocity_npz (e.g. 'v_lift_R5'): per-galaxy velocities shared by Liam Parker
        (ksz_lrg_dr2_20261117, ksz_bgs_dr2_20261119). Their targetid column was rounded through float64, so
        objects are matched on the exact TARGETID read from the clustering catalogues at source_row.
    any column of velocity_fits (e.g. 'vR'), when source.velocity_fits is set instead of velocity_npz:
        per-cap FITS catalogues without TARGETID (Boryana Hadzhiyska's reconstructions), matched on the exact
        (RA, DEC, Z) of each object.

Run from desi/spec_Y3/catalogue/ in the cosmodesi environment:
    python -u swap_velocities.py -p ./configs/velocity_dr2/prepare_cat_LRG_zeroV_v_lift_R5.yaml
"""
import argparse
import os

import numpy as np
import pandas as pd
import yaml
from astropy.table import Table
from cosmoprimo.fiducial import DESI # type: ignore


def exact_velocity_table(velocity_npz: str, clustering_fits: str, column: str) -> pd.DataFrame:
    """Per-galaxy velocities indexed by the exact TARGETID.

    Args:
        velocity_npz (str): galaxies.npz with cap, source_row, targetid, z and the velocity column.
        clustering_fits (str): clustering catalogue path with a '{cap}' placeholder (NGC, SGC).
        column (str): velocity column to return.

    Returns:
        pandas.DataFrame: columns Z_VEL and VEL, indexed by TARGETID.
    """
    g = np.load(velocity_npz)
    cap, source_row = g['cap'], g['source_row']
    tid = np.zeros(len(cap), dtype=np.int64)
    for c, name in enumerate(['NGC', 'SGC']):
        sel = cap == c
        t = Table.read(clustering_fits.format(cap=name), format='fits')
        tid[sel] = np.asarray(t['TARGETID'])[source_row[sel]]
        if not np.array_equal(np.asarray(t['Z'])[source_row[sel]], g['z'][sel]):
            raise ValueError(f"{name}: redshifts at source_row differ from {velocity_npz}")
    # the shared targetid is the exact one rounded through float64
    if not np.array_equal(tid.astype(np.float64).astype(np.int64), g['targetid']):
        raise ValueError("TARGETID at source_row does not round to the shared targetid")
    if len(np.unique(tid)) != len(tid):
        raise ValueError("TARGETID not unique in the velocity file")
    return pd.DataFrame({'Z_VEL': g['z'], 'VEL': g[column].astype(np.float64)}, index=pd.Index(tid, name='TARGETID'))


def coordinate_velocity_table(velocity_fits: str, column: str) -> pd.DataFrame:
    """Per-galaxy velocities indexed by the exact (RA, DEC, Z), for catalogues without TARGETID.

    Args:
        velocity_fits (str): catalogue path with a '{cap}' placeholder (NGC, SGC), with RA, DEC, Z and the
            velocity column.
        column (str): velocity column to return.

    Returns:
        pandas.DataFrame: columns Z_VEL and VEL, indexed by (RA, DEC, Z).
    """
    frames = []
    for name in ['NGC', 'SGC']:
        t = Table.read(velocity_fits.format(cap=name), format='fits')
        frames.append(pd.DataFrame({key: np.asarray(t[key], dtype=np.float64) for key in ['RA', 'DEC', 'Z', column]}))
    v = pd.concat(frames, ignore_index=True)
    v = pd.DataFrame({'Z_VEL': v['Z'].values, 'VEL': v[column].values},
                     index=pd.MultiIndex.from_arrays([v['RA'].values, v['DEC'].values, v['Z'].values]))
    if not v.index.is_unique:
        raise ValueError(f"(RA, DEC, Z) not unique in {velocity_fits}")
    return v


parser = argparse.ArgumentParser(description='Swap the velocities of the production catalogues.')
parser.add_argument('-p', '--path2config', type=str, required=True, help='Path to the configuration file.')
parser.add_argument('--overwrite', action='store_true', help='Overwrite existing output catalogues.')
args = vars(parser.parse_args())
print(f"Arguments: {args}")
with open(args['path2config'], 'r') as f:
    config = yaml.safe_load(f)

production_dir = config['source']['production_dir']
column = config['source']['velocity_column']
filter_type = config['processing'].get('filter_type', 'no_src_with_cluster_mask')
renorm = config['processing'].get('renorm', False)
z_eff = config['processing'].get('z_eff', 0.780)
zeroV = config['processing'].get('zeroV', True)
zbins = [tuple(z) for z in config['processing']['zbins']]
fields = config['processing'].get('fields', ['full'])
output_dir = config['save']['save_dir']

if os.path.realpath(output_dir) == os.path.realpath(production_dir):
    raise ValueError("save_dir must differ from production_dir")
os.makedirs(output_dir, exist_ok=True)
os.chmod(output_dir, 0o700) # unpublished velocities (LiFT, collaborators' reconstructions): private

if 'velocity_fits' in config['source'] and 'velocity_npz' in config['source']:
    raise ValueError("set source.velocity_fits or source.velocity_npz, not both")
match_coordinates = 'velocity_fits' in config['source']
if column != 'official':
    if match_coordinates:
        velocities = coordinate_velocity_table(config['source']['velocity_fits'], column)
        print(f"{column}: {len(velocities)} galaxies with exact (RA, DEC, Z)")
    else:
        velocities = exact_velocity_table(config['source']['velocity_npz'], config['source']['clustering_fits'], column)
        print(f"{column}: {len(velocities)} galaxies with exact TARGETID")
cosmo = DESI()

names = [f"{field}_zbin{i}_{filter_type}.csv" for field in fields for i in range(1, len(zbins) + 1)]
existing = [n for n in names if os.path.exists(os.path.join(output_dir, n))]
if existing and not args['overwrite']:
    raise FileExistsError(f"{existing} exist in {output_dir} (use --overwrite)")

for field in fields:
    for i in range(1, len(zbins) + 1):
        name = f"{field}_zbin{i}_{filter_type}.csv"
        output_file = os.path.join(output_dir, name)
        # round_trip: the other columns are copied bit for bit
        df = pd.read_csv(os.path.join(production_dir, name), float_precision='round_trip')

        if column == 'official':
            vel = df['VEL_LOS'].values.astype(np.float64)
        else:
            if match_coordinates:
                matched = velocities.reindex(pd.MultiIndex.from_arrays([df['RA'].values, df['DEC'].values, df['Z'].values]))
            else:
                matched = velocities.reindex(df['TARGETID'].values)
            if matched['VEL'].isna().any():
                raise ValueError(f"{name}: {matched['VEL'].isna().sum()} objects without a {column} velocity")
            # (always 0 on the coordinate route, where Z is part of the key)
            dz = np.max(np.abs(matched['Z_VEL'].values - df['Z'].values))
            if dz > 1e-10:
                raise ValueError(f"{name}: redshifts differ from the velocity file by up to {dz:.2e}")
            vel = matched['VEL'].values
            df = df.drop(columns=['DISP_LOS']) # the official displacement, not this velocity's

        vel_renorm = vel * cosmo.growth_factor(df['Z'].values) / cosmo.growth_factor(z_eff) if renorm else vel.copy()
        # zeroV: subtract the mean of each cap within the bin, as prepare_catalogue_yaml.py does
        if zeroV:
            for gc in np.unique(df['GC']):
                sel = (df['GC'] == gc).values
                vel_renorm[sel] -= np.mean(vel_renorm[sel])

        if column == 'official': # known answer, checked on every run
            np.testing.assert_allclose(vel_renorm, df['VEL_LOS_RENORM'].values, rtol=0, atol=1e-8,
                                       err_msg=f"{name}: official velocities do not reproduce the production VEL_LOS_RENORM")
        df['VEL_LOS'] = vel
        df['VEL_LOS_RENORM'] = vel_renorm
        # temporary file first, so a failed write never leaves a truncated catalogue
        df.to_csv(output_file + '.tmp', index=False)
        os.replace(output_file + '.tmp', output_file)
        print(f"Saved {name}: {len(df)} objects, VEL_LOS_RENORM std {np.std(vel_renorm):.1f} km/s -> {output_file}")

print('done!!')
