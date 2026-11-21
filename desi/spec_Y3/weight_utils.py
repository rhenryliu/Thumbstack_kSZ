"""
weight_utils.py
Per-object completeness weights for the kSZ stacks (ThumbStack estimator 'ksz_compweight').

Used by make_stacks_yaml.py and make_stacks_shuffle.py. The weights are the product of catalogue
columns (default WEIGHT_COMP * WEIGHT_ZFAIL), read from a table that is row aligned with the
ThumbStack catalogue: the intermediate CSV next to it (as for the mass filter), or a side file
(e.g. the BGS CIGALE lens weights). Alternatively they are joined by TARGETID from other catalogues
(e.g. the DESI v2 PIP clustering catalogues); objects not found get NaN, which ThumbStack's
catalogue mask drops.

load_mstar_weights gives the stellar-mass weights of the estimator 'ksz_mstarweight', w = M*^power,
from the log10 stellar mass column of the same row-aligned table.
"""

from fractions import Fraction

import numpy as np
import pandas as pd

from mass_filter_utils import get_csv_path

DEFAULT_COLUMNS = ['WEIGHT_COMP', 'WEIGHT_ZFAIL']
# RA, DEC [deg] and Z of the weight table must equal the catalogue's to this tolerance
POSITION_TOLERANCE = 1.e-8


def read_table(path: str, columns: list) -> pd.DataFrame:
    """Read *columns* of a CSV or FITS table.

    Args:
        path: Path to a ``.csv`` or ``.fits`` file.
        columns: Column names to read.

    Returns:
        DataFrame with those columns, in file order.
    """
    if path.endswith('.csv'):
        return pd.read_csv(path, usecols=columns)[columns]
    import fitsio
    data = fitsio.read(path, columns=columns)
    return pd.DataFrame({c: data[c].astype(np.float64) if c != 'TARGETID' else data[c] for c in columns})


def check_alignment(table: pd.DataFrame, catalog, path: str) -> None:
    """Check that *table* is row aligned with the ThumbStack catalogue.

    Args:
        table: Weight table with RA, DEC and Z columns.
        catalog: ThumbStack Catalog (RA, DEC, Z, nObj).
        path: Path of the table, for the error message.

    Raises:
        ValueError: If the lengths differ or any RA, DEC or Z differs by more than
            POSITION_TOLERANCE.
    """
    if len(table) != catalog.nObj:
        raise ValueError(f"[weights] {len(table)} rows in {path}, {catalog.nObj} objects in the catalogue")
    for col, ref in (('RA', catalog.RA), ('DEC', catalog.DEC), ('Z', catalog.Z)):
        diff = np.max(np.abs(table[col].values.astype(np.float64) - ref))
        if not diff <= POSITION_TOLERANCE:
            raise ValueError(f"[weights] {col} of {path} differs from the catalogue by up to {diff}")


def join_by_targetid(targetid: np.ndarray, paths: list, columns: list) -> np.ndarray:
    """Product of *columns* from the catalogues in *paths*, matched by TARGETID.

    Args:
        targetid: TARGETID of each catalogue object.
        paths: Catalogues holding TARGETID and *columns* (e.g. the PIP NGC and SGC files).
        columns: Weight columns to multiply.

    Returns:
        One weight per object; NaN where the TARGETID is not found.

    Raises:
        ValueError: If a TARGETID appears more than once in the catalogues of *paths*.
    """
    # int64 throughout: mixing int64 and uint64 would promote to float64, where neighbouring
    # TARGETIDs (~4e16 > 2^53) compare equal
    tables = [read_table(p, ['TARGETID'] + columns) for p in paths]
    for tab, p in zip(tables, paths):
        if not np.issubdtype(tab['TARGETID'].dtype, np.integer):
            raise ValueError(f"[weights] TARGETID of {p} is {tab['TARGETID'].dtype}, not integer")
        tab['TARGETID'] = tab['TARGETID'].values.astype(np.int64)
    if not np.issubdtype(np.asarray(targetid).dtype, np.integer):
        raise ValueError(f"[weights] TARGETID of the catalogue is {np.asarray(targetid).dtype}, not integer")
    targetid = np.asarray(targetid).astype(np.int64)
    other = pd.concat(tables, ignore_index=True)
    if other['TARGETID'].duplicated().any():
        raise ValueError(f"[weights] duplicate TARGETID in {paths}")
    w_other = np.prod(other[columns].values.astype(np.float64), axis=1)
    order = np.argsort(other['TARGETID'].values)
    tid_sorted = other['TARGETID'].values[order]
    i = np.clip(np.searchsorted(tid_sorted, targetid), 0, len(tid_sorted) - 1)
    found = tid_sorted[i] == targetid
    w = np.full(len(targetid), np.nan)
    w[found] = w_other[order][i][found]
    return w


def load_weights(weights_config: dict, catalog, cat_name: str, cat_type: str, cat_dir: str,
                 field: str, filter_cut: str) -> np.ndarray:
    """Completeness weight of each catalogue object, as configured in ``stack.weights``.

    Args:
        weights_config: ``columns`` (list, default WEIGHT_COMP and WEIGHT_ZFAIL); ``file``
            (optional table row aligned with the catalogue, default the intermediate CSV given by
            get_csv_path); ``join`` (optional list of catalogues to take the columns from, matched
            by the TARGETID of the row-aligned table).
        catalog: ThumbStack Catalog the weights are for.
        cat_name, cat_type, cat_dir, field, filter_cut: As for get_csv_path.

    Returns:
        Array of length ``catalog.nObj``; NaN for objects without a weight (``join`` only).

    Raises:
        ValueError: If the table is not row aligned with the catalogue, or a weight is not finite
            and positive where one is expected.
    """
    columns = list(weights_config.get('columns', DEFAULT_COLUMNS))
    path = weights_config.get('file') or get_csv_path(cat_name, cat_type, cat_dir, field, filter_cut)
    join = weights_config.get('join')
    table = read_table(path, ['RA', 'DEC', 'Z'] + (['TARGETID'] if join else columns))
    check_alignment(table, catalog, path)
    if join:
        w = join_by_targetid(table['TARGETID'].values, list(join), columns)
        source = f"{columns} joined by TARGETID from {list(join)}"
    else:
        w = np.prod(table[columns].values.astype(np.float64), axis=1)
        source = f"{columns} of {path}"
    found = np.isfinite(w)
    if not found.any():
        raise ValueError(f"[weights] no object has a weight: {source}")
    if np.any(w[found] <= 0.) or (not join and not np.all(found)):
        raise ValueError(f"[weights] weights not finite and positive: {source}")
    print(f"[weights] {cat_name}: {source}; {found.sum()} of {len(w)} objects with a weight, "
          f"mean {np.mean(w[found]):.4f}, N_eff/N {np.sum(w[found])**2 / np.sum(w[found]**2) / found.sum():.4f}")
    return w


def load_mstar_weights(mstar_config: dict, catalog, cat_name: str, cat_type: str, cat_dir: str,
                       field: str, filter_cut: str) -> np.ndarray:
    """Stellar-mass weight of each catalogue object, w = (M*)^power, or with a floor
    w = (M*)^power + (M_min)^power, as configured in ``stack.mstar_weights`` or in one entry of
    ``stack.mstar_weightings``.

    M* is the linear mass, 10^col. The normalisation of w cancels in the weighted estimator, so the
    mass units (and factors of h) do not matter, as long as the floor is in the units of col. Objects with log10 M* <= 0 (failed fits) or not finite get NaN, which ThumbStack's
    catalogue mask drops from every estimator of the run.

    Args:
        mstar_config: ``power`` (required; a number or a fraction string such as ``'2/3'``);
            ``col`` (log10 stellar mass column, default LOGMSTAR); ``file`` (optional table row
            aligned with the catalogue, default the intermediate CSV given by get_csv_path);
            ``floor_logmstar`` (optional log10 M_min: a number, or a dict keyed by catalogue
            name; without it, no floor).
        catalog: ThumbStack Catalog the weights are for.
        cat_name, cat_type, cat_dir, field, filter_cut: As for get_csv_path.

    Returns:
        Array of length ``catalog.nObj``; NaN for objects without a valid stellar mass.

    Raises:
        ValueError: If the table is not row aligned with the catalogue, no object has a valid
            stellar mass, or ``floor_logmstar`` has no value for *cat_name*.
    """
    power = float(Fraction(str(mstar_config['power'])))
    floor = mstar_config.get('floor_logmstar')
    if isinstance(floor, dict):
        if cat_name not in floor:
            raise ValueError(f"[mstar weights] floor_logmstar has no value for {cat_name}")
        floor = floor[cat_name]
    col = mstar_config.get('col', 'LOGMSTAR')
    path = mstar_config.get('file') or get_csv_path(cat_name, cat_type, cat_dir, field, filter_cut)
    table = read_table(path, ['RA', 'DEC', 'Z', col])
    check_alignment(table, catalog, path)
    logm = table[col].values.astype(np.float64)
    found = np.isfinite(logm) & (logm > 0.)
    if not found.any():
        raise ValueError(f"[mstar weights] no object has a valid {col} in {path}")
    w = np.full(len(logm), np.nan)
    w[found] = 10.**(power * logm[found])
    if floor is not None:
        w[found] += 10.**(power * float(floor))
    floor_str = '' if floor is None else f" + M_min^{power:.6g}, log10 M_min {float(floor):.4f}"
    print(f"[mstar weights] {cat_name}: M*^{power:.6g}{floor_str} from {col} of {path}; {found.sum()} of {len(w)} objects "
          f"with {col} > 0, N_eff/N {np.sum(w[found])**2 / np.sum(w[found]**2) / found.sum():.4f}")
    return w
