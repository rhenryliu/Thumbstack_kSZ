"""
mass_filter_utils.py
Utilities for applying stellar mass filters to catalogues at analysis time.

Used by make_stacks_yaml.py on the kSZ side. A partial copy (compute_mass_mask
and mass_label only) lives in DESIxHSC-Lensing/src/mass_filter_utils.py for
the lensing precompute pipeline.

Not sample-specific: works for any catalogue that has a stellar mass column in
the intermediate CSV, regardless of whether it is BGS, LRG, or another type.
"""

import os
import numpy as np


def compute_mass_mask(logmstar: np.ndarray, mass_config: dict) -> np.ndarray:
    """Return a boolean mask selecting objects by stellar mass.

    Args:
        logmstar: 1-D array of log10 stellar masses. NaN values are always
            excluded from any selection.
        mass_config: Dict controlling the selection strategy. Required key:
            ``strategy`` (str). Additional keys depend on the strategy:

            ``'top_frac'``
                Keep the top *fraction* of objects by mass.
                Required: ``fraction`` (float, e.g. 0.05 for top 5%).

            ``'percentile_bins'``
                Divide the sample into *n_bins* equal-count bins and select
                the bin at *bin_idx* (0 = lowest mass, n_bins-1 = highest).
                Required: ``n_bins`` (int), ``bin_idx`` (int).

            ``'cumulative'``
                Keep all objects with LOGMSTAR >= *threshold*.
                Required: ``threshold`` (float).

            ``'range'``
                Keep objects with *min* <= LOGMSTAR < *max*.
                Required: ``min`` (float), ``max`` (float).

    Returns:
        Boolean array of the same length as *logmstar*.

    Raises:
        ValueError: If *strategy* is not one of the recognised values.
    """
    logmstar = np.asarray(logmstar, dtype=float)
    strategy = mass_config['strategy']
    finite = np.isfinite(logmstar)

    if strategy == 'top_frac':
        fraction = float(mass_config['fraction'])
        cutoff = np.nanpercentile(logmstar, (1.0 - fraction) * 100.0)
        mask = finite & (logmstar >= cutoff)

    elif strategy == 'percentile_bins':
        n_bins = int(mass_config['n_bins'])
        bin_idx = int(mass_config['bin_idx'])
        edges = np.nanpercentile(logmstar, np.linspace(0, 100, n_bins + 1))
        lo = edges[bin_idx]
        hi = edges[bin_idx + 1]
        if bin_idx == n_bins - 1:
            # Include the maximum-mass object in the top bin.
            mask = finite & (logmstar >= lo)
        else:
            mask = finite & (logmstar >= lo) & (logmstar < hi)

    elif strategy == 'cumulative':
        threshold = float(mass_config['threshold'])
        mask = finite & (logmstar >= threshold)

    elif strategy == 'range':
        lo = float(mass_config['min'])
        hi = float(mass_config['max'])
        mask = finite & (logmstar >= lo) & (logmstar < hi)

    else:
        raise ValueError(
            f"Unknown mass filter strategy: {strategy!r}. "
            "Must be one of: 'top_frac', 'percentile_bins', 'cumulative', 'range'."
        )

    selected = logmstar[mask]
    log_mean = np.log10(np.mean(10.0 ** selected))
    print(f"[mass filter] log10(M*): min={selected.min():.3f}, "
          f"mean={log_mean:.3f}, max={selected.max():.3f}")
    return mask


def mass_label(mass_config: dict) -> str:
    """Return a short, filesystem-safe label string for a mass filter config.

    This is the single source of truth for naming on both the kSZ and lensing
    sides. Both pipelines call this function so that path suffixes are always
    consistent.

    Args:
        mass_config: Same dict as passed to :func:`compute_mass_mask`.

    Returns:
        A concise label string. Examples::

            top_frac  fraction=0.05          ->  'top5pct'
            top_frac  fraction=0.01          ->  'top1pct'
            percentile_bins n=5 idx=4        ->  'mbin4of5'
            percentile_bins n=5 idx=0        ->  'mbin0of5'
            cumulative  threshold=10.5       ->  'mge10p5'
            cumulative  threshold=11.0       ->  'mge11p0'
            range  min=10.0 max=10.5         ->  'm10p0to10p5'

    Raises:
        ValueError: If *strategy* is not one of the recognised values.
    """
    strategy = mass_config['strategy']

    if strategy == 'top_frac':
        pct = int(round(float(mass_config['fraction']) * 100))
        return f'top{pct}pct'

    elif strategy == 'percentile_bins':
        n = int(mass_config['n_bins'])
        idx = int(mass_config['bin_idx'])
        return f'mbin{idx}of{n}'

    elif strategy == 'cumulative':
        val_str = f"{float(mass_config['threshold']):.1f}".replace('.', 'p')
        return f'mge{val_str}'

    elif strategy == 'range':
        lo = f"{float(mass_config['min']):.1f}".replace('.', 'p')
        hi = f"{float(mass_config['max']):.1f}".replace('.', 'p')
        return f'm{lo}to{hi}'

    else:
        raise ValueError(
            f"Unknown mass filter strategy: {strategy!r}. "
            "Must be one of: 'top_frac', 'percentile_bins', 'cumulative', 'range'."
        )


def get_csv_path(cat_name: str, cat_type: str, cat_dir: str,
                 field: str, filter_cut: str) -> str:
    """Derive the intermediate CSV path for a ThumbStack catalog name.

    The intermediate CSVs produced by prepare_catalogue_all_yaml.py or
    prepare_catalogue_yaml.py live in *cat_dir* alongside the ThumbStack
    catalog subdirectories. This function encodes the fixed naming convention::

        DESIY3_{cat_type}       ->  {cat_dir}/{field}_catalog_Y3_{filter_cut}.csv
        DESIY3_{cat_type}_z{n}  ->  {cat_dir}/{field}_zbin{n}_{filter_cut}.csv

    Works for any cat_type string (BGS, LRG, etc.).

    Args:
        cat_name: ThumbStack catalog name, e.g. ``'DESIY3_BGS'`` or
            ``'DESIY3_LRG_z3'``.
        cat_type: Sample type string, e.g. ``'BGS'`` or ``'LRG'``.
        cat_dir: Directory containing the intermediate CSVs and ThumbStack
            subdirectories (value of ``data.cat_dir`` in the stacking YAML).
        field: Galactic field: ``'full'``, ``'NGC'``, or ``'SGC'``.
        filter_cut: Filter type, e.g. ``'nopairs'``.

    Returns:
        Absolute path to the intermediate CSV file.

    Raises:
        ValueError: If *cat_name* does not match the expected pattern.
    """
    base = f'DESIY3_{cat_type}'

    if cat_name == base:
        csv_name = f'{field}_catalog_Y3_{filter_cut}.csv'
    elif cat_name.startswith(base + '_z'):
        zbin = cat_name[len(base) + 2:]  # e.g. '1', '2', '3'
        csv_name = f'{field}_zbin{zbin}_{filter_cut}.csv'
    else:
        raise ValueError(
            f"Cannot derive CSV path for catalog name {cat_name!r}. "
            f"Expected 'DESIY3_{cat_type}' or 'DESIY3_{cat_type}_z{{n}}'."
        )

    return os.path.join(cat_dir, csv_name)


def apply_catalog_mask(catalog, mask: np.ndarray) -> None:
    """Apply a boolean mask in-place to all arrays in a ThumbStack Catalog.

    Applies *mask* to every array attribute written by Catalog.writeCatalog()
    and read by Catalog.loadCatalog() in catalog.py (24 columns total). Fields
    absent from the object are silently skipped. Updates catalog.nObj.

    Args:
        catalog: ThumbStack Catalog object. Modified in-place.
        mask: Boolean array of length ``catalog.nObj``.
    """
    _ARRAY_FIELDS = [
        'RA', 'DEC', 'Z',
        'coordX', 'coordY', 'coordZ',
        'dX', 'dY', 'dZ',
        'dXKaiser', 'dYKaiser', 'dZKaiser',
        'vX', 'vY', 'vZ',
        'vR', 'vTheta', 'vPhi',
        'Mstellar', 'hasM', 'Mvir',
        'integratedTau', 'integratedKSZ', 'integratedY',
    ]
    mask = np.asarray(mask, dtype=bool)
    for attr in _ARRAY_FIELDS:
        if hasattr(catalog, attr):
            setattr(catalog, attr, getattr(catalog, attr)[mask])
    catalog.nObj = int(mask.sum())
