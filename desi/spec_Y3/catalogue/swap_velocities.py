"""Swap the velocities of the production LRG (or BGS) catalogues for another velocity estimate.

The production per-bin catalogues ({field}_zbin{i}_{filter_type}.csv, written by prepare_catalogue_yaml.py)
are kept object for object and in the same row order, and only VEL_LOS and VEL_LOS_RENORM are replaced, so
the stacks can reuse the production per-object filtering (stack.filtering_dir in make_stacks_yaml.py).
convert_to_TS_catalog.py then writes the ThumbStack catalogues from the same config.

velocity_column:
    'official': the production VEL_LOS (BAO RecSym reconstruction), the known-answer control: with
        renorm and zeroV it must reproduce the production VEL_LOS_RENORM (re-zeroed per zero_within group for
        the bins listed there).

processing.zero_within (optional, {bin: [[z_lo, z_hi], ...]} tiling that bin): zero the bin per cap within each
sub-range z_lo < z <= z_hi instead of over the whole bin.
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


def zero_groups(v: np.ndarray, groups: np.ndarray) -> np.ndarray:
    """Copy of v with the mean of each group subtracted.

    Args:
        v (numpy.ndarray): velocities.
        groups (numpy.ndarray): group label of each object (e.g. its cap).

    Returns:
        numpy.ndarray: the zeroed velocities.
    """
    v = v.copy()
    for g in np.unique(groups):
        sel = groups == g
        v[sel] -= np.mean(v[sel])
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
# optional {bin index: [[z_lo, z_hi], ...]}: zero that bin per cap within each sub-range (z_lo < z <= z_hi)
zero_within = {int(b): ranges for b, ranges in config['processing'].get('zero_within', {}).items()}
if zero_within and not zeroV:
    raise ValueError("processing.zero_within needs zeroV: true")
for b, ranges in zero_within.items():
    if b not in range(1, len(zbins) + 1):
        raise ValueError(f"zero_within bin {b} is not in 1..{len(zbins)}")
    if any(hi != lo for (_, hi), (lo, _) in zip(ranges[:-1], ranges[1:])) or (ranges[0][0], ranges[-1][1]) != zbins[b - 1]:
        raise ValueError(f"zero_within[{b}] = {ranges} must tile bin {b} = {zbins[b - 1]} without gaps or overlaps")
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
        # zeroV: subtract the mean of each cap within the bin, as prepare_catalogue_yaml.py does;
        # with zero_within, of each cap within each listed sub-range of the bin
        groups = df['GC'].values.astype(str)
        if i in zero_within:
            sub = np.full(len(df), -1)
            for k, (lo, hi) in enumerate(zero_within[i]):
                sub[(df['Z'].values > lo) & (df['Z'].values <= hi)] = k
            if (sub < 0).any():
                raise ValueError(f"{name}: {np.sum(sub < 0)} objects outside the zero_within sub-ranges")
            groups = np.char.add(np.char.add(groups, '|'), sub.astype(str))
        if zeroV:
            vel_renorm = zero_groups(vel_renorm, groups)
            if max(abs(np.mean(vel_renorm[groups == g])) for g in np.unique(groups)) > 1e-8:
                raise ValueError(f"{name}: a zeroing group has a non-zero mean velocity")

        if column == 'official': # known answer, checked on every run
            # the production VEL_LOS_RENORM is zeroed per cap over the bin; re-zero it per sub-range to compare
            reference = zero_groups(df['VEL_LOS_RENORM'].values, groups) if i in zero_within else df['VEL_LOS_RENORM'].values
            np.testing.assert_allclose(vel_renorm, reference, rtol=0, atol=1e-8,
                                       err_msg=f"{name}: official velocities do not reproduce the production VEL_LOS_RENORM")
        df['VEL_LOS'] = vel
        df['VEL_LOS_RENORM'] = vel_renorm
        # temporary file first, so a failed write never leaves a truncated catalogue
        df.to_csv(output_file + '.tmp', index=False)
        os.replace(output_file + '.tmp', output_file)
        print(f"Saved {name}: {len(df)} objects, VEL_LOS_RENORM std {np.std(vel_renorm):.1f} km/s -> {output_file}")

print('done!!')
