import numpy as np
from astropy.coordinates import SkyCoord
import astropy.units as u
from astropy.table import Table
from astropy.cosmology import Planck18 

import sys
sys.path.append('../src/')
# from utils import arcmin_to_rp, rp_to_arcmin

def rp_to_arcmin(rp_mpc_h, z_eff, cosmology=Planck18, comoving=True):
    """Convert transverse separation rp [Mpc/h] into angular scale [arcmin].

    Args:
        rp_mpc_h (float or array-like): Transverse separation in Mpc/h. 
            Interpreted as comoving if `comoving=True`, physical if `comoving=False`.
        z_eff (float): Effective redshift.
        cosmology (astropy.cosmology, optional): Cosmology object. Defaults to Planck18.
        comoving (bool, optional): If True, use transverse comoving distance D_M(z).
            If False, use angular diameter distance D_A(z). Defaults to False.

    Returns:
        float or array-like: Angular scale in arcminutes.
    """

    h = cosmology.h
    rp_mpc = np.asarray(rp_mpc_h) / h  # convert Mpc/h → Mpc

    if comoving:
        dist = cosmology.comoving_transverse_distance(z_eff)  # D_M(z)
    else:
        dist = cosmology.angular_diameter_distance(z_eff)     # D_A(z)

    theta = (rp_mpc * u.Mpc) / dist * u.rad
    return theta.to(u.arcmin).value

def arcmin_to_rp(theta_arcmin, z_eff, cosmology=Planck18, comoving=True):
    """Convert angular scale [arcmin] into transverse separation rp [Mpc/h].

    Args:
        theta_arcmin (float or array-like): Angular separation in arcminutes.
        z_eff (float): Effective redshift.
        cosmology (astropy.cosmology, optional): Cosmology object. Defaults to Planck18.
        comoving (bool, optional): If True, use transverse comoving distance D_M(z).
            If False, use angular diameter distance D_A(z). Defaults to False.

    Returns:
        float or array-like: Transverse separation in Mpc/h. 
            Comoving if `comoving=True`, physical if `comoving=False`.
    """
    h = cosmology.h
    theta_rad = np.asarray(theta_arcmin) * (u.arcmin).to(u.rad)  # arcmin → rad

    if comoving:
        dist = cosmology.comoving_transverse_distance(z_eff)  # D_M(z)
    else:
        dist = cosmology.angular_diameter_distance(z_eff)     # D_A(z)

    rp_mpc = (theta_rad * dist).to(u.Mpc).value  # Mpc
    rp_mpc_h = rp_mpc * h  # Mpc/h
    return rp_mpc_h

def enforce_min_separation(
    tab: Table,
    ra_col: str = "RA",
    dec_col: str = "Dec",
    min_sep = 1.0 * u.arcsec,
    frame: str = "icrs",
    priority_col = None,
    prefer_higher_priority: bool = True,
    return_mask: bool = False,
) :
    """
    Return a sub-table where no two rows are within `min_sep` on the sky.

    Parameters
    ----------
    tab : astropy.table.Table
        Input catalog containing sky positions (RA/Dec).
    ra_col, dec_col : str
        Column names for right ascension and declination (degrees by default).
    min_sep : astropy.units.Quantity or float
        Minimum allowed separation. If Quantity, must have angle units.
        If float (unitless), interpreted as arcsec.
    frame : str
        Astropy frame string for the coordinates (e.g. "icrs", "fk5").
    priority_col : str or None
        Optional column name giving each row's priority. The greedy solver
        keeps higher (or lower) values first depending on `prefer_higher_priority`.
        If None, keeps rows in their original order.
    prefer_higher_priority : bool
        If True, larger values in `priority_col` are preferred. If False, smaller.
    return_mask : bool
        If True, also return a boolean mask of kept rows with length len(tab).

    Returns
    -------
    sub : astropy.table.Table
        Pruned catalog with min separation enforced.
    keep_mask : np.ndarray (optional)
        Boolean mask indexing the kept rows from the original table.

    Notes
    -----
    - This uses a greedy maximal independent set: it guarantees no pair is
      closer than `min_sep` but does not guarantee the *maximum* possible count.
      In practice, with a sensible priority (e.g., brightness), it’s what you want.
    - Rows with NaN or masked RA/Dec are *kept* as-is (they don't participate in
      the distance checks); change behavior below if you want to drop them.
    """
    if ra_col not in tab.colnames or dec_col not in tab.colnames:
        raise KeyError(f"Table must contain '{ra_col}' and '{dec_col}' columns.")

    # Normalize min_sep
    if not isinstance(min_sep, u.Quantity):
        min_sep = float(min_sep) * u.arcsec
    if not min_sep.unit.is_equivalent(u.deg): # type: ignore
        raise ValueError("min_sep must be an angular quantity (e.g. arcsec, arcmin, deg).")

    # Identify valid (finite) coordinates; we’ll allow invalid rows to pass through
    ra = np.array(tab[ra_col], dtype=float)
    dec = np.array(tab[dec_col], dtype=float)
    valid = np.isfinite(ra) & np.isfinite(dec)

    N = len(tab)
    keep_mask = np.zeros(N, dtype=bool)

    # Build coordinates for valid rows
    if valid.any():
        coords = SkyCoord(ra=ra[valid] * u.deg, dec=dec[valid] * u.deg, frame=frame)

        # Find all pairs within min_sep (including self-matches)
        i, j, sep2d, _ = coords.search_around_sky(coords, min_sep)

        # Remove self-matches and duplicate orderings; keep i < j
        pair_mask = i < j
        i = i[pair_mask]
        j = j[pair_mask]

        # Map back to original indices
        valid_idx = np.nonzero(valid)[0]
        i_full = valid_idx[i]
        j_full = valid_idx[j]

        # Build adjacency list (neighbors within min_sep) for greedy selection
        neighbors = [[] for _ in range(N)]
        for a, b in zip(i_full, j_full):
            neighbors[a].append(b)
            neighbors[b].append(a)

        # Order to process: by priority if provided, else original index
        if priority_col is not None:
            if priority_col not in tab.colnames:
                raise KeyError(f"priority_col '{priority_col}' not found in table.")
            pr = np.array(tab[priority_col])
            # Replace NaN with -inf/+inf depending on preference so they sort last
            pr = pr.astype(float, copy=False)
            nan_fill = -np.inf if prefer_higher_priority else np.inf
            pr = np.where(np.isfinite(pr), pr, nan_fill)
            order = np.argsort(pr)  # ascending
            if prefer_higher_priority:
                order = order[::-1]  # descending
        else:
            order = np.arange(N)

        blocked = np.zeros(N, dtype=bool)
        for idx in order:
            if blocked[idx]:
                continue
            # Accept this source
            keep_mask[idx] = True
            # Block all neighbors from being selected later
            for nb in neighbors[idx]:
                blocked[nb] = True

    else:
        # No valid coordinates; keep everything (or change to keep nothing)
        keep_mask[:] = True

    # Ensure rows with invalid positions are kept (as documented)
    keep_mask |= ~valid

    sub = tab[keep_mask]
    return (sub, keep_mask) if return_mask else sub # type: ignore



def mean_redshift_from_table(
    tab: Table,
    z_col: str = "Z",
    weights_col = None,
) -> float:
    """
    Compute a robust mean redshift from a catalog.

    - Ignores non-finite z.
    - Optional weights (e.g., completeness or inverse-variance) supported.
    """
    if z_col not in tab.colnames:
        raise KeyError(f"Table must contain a redshift column '{z_col}'.")
    z = np.asarray(tab[z_col], dtype=float)
    m = np.isfinite(z)
    if not m.any():
        raise ValueError("No finite redshift values found.")
    if weights_col is None:
        return float(np.mean(z[m]))
    if weights_col not in tab.colnames:
        raise KeyError(f"weights_col '{weights_col}' not found in table.")
    w = np.asarray(tab[weights_col], dtype=float)
    w = np.where(np.isfinite(w), w, 0.0)
    w = w[m]
    z = z[m]
    if np.sum(w) <= 0:
        raise ValueError("All weights are zero or non-finite.")
    return float(np.average(z, weights=w))

def physical_to_min_sep(
    d_target,
    z_mean: float,
    *,
    cosmology=Planck18,
    comoving: bool = False,
    out_unit: u.Unit = u.arcsec,
):
    """
    Convert a target transverse separation into an angular min_sep at z_mean.

    Parameters
    ----------
    d_target : astropy.units.Quantity
        Target transverse separation. If comoving=False, this should be a *physical*
        (proper) distance (e.g., 50 * u.kpc). If comoving=True, pass a *comoving*
        distance (e.g., 50 * u.kpc) and it will be divided by (1+z_mean).
    z_mean : float
        Mean redshift at which to enforce the separation.
    cosmology : astropy.cosmology.Cosmology
        Cosmology to use (default: Planck18).
    comoving : bool
        If True, interpret `d_target` as comoving and convert to physical.
    out_unit : astropy.units.Unit
        Angular unit for the returned angle (default: arcsec).

    Returns
    -------
    theta_min : astropy.units.Quantity
        Angular separation corresponding to `d_target` at `z_mean`.
    """
    if not isinstance(d_target, u.Quantity):
        raise TypeError("d_target must be an astropy Quantity with length units.")
    if not d_target.unit.is_equivalent(u.kpc): # type: ignore
        raise ValueError("d_target must have length units (e.g., kpc, Mpc).")

    # physical (proper) distance at z_mean
    d_phys = d_target.to(u.Mpc)
    if comoving:
        d_phys = d_phys / (1.0 + z_mean)

    # angular diameter distance (proper per radian)
    DA = cosmology.angular_diameter_distance(z_mean)  # in Mpc
    # theta = (d_phys / DA).to(u.rad)  # exact small-angle geometry
    theta = (d_phys / DA) * (u.rad)  # exact small-angle geometry

    return theta.to(out_unit)
    

if __name__ == "__main__":
    h = Planck18.h
    
    redshifts = [0.47461575, 0.6352965, 0.79460454, 0.92029834] # mean redshifts for photo-z bins 1-4
    physical_sep = 1.0 * u.Mpc# / h  # target minimum physical (proper distance)
    
    # bin_num = 1
    bins = [1, 2, 3, 4]
    for bin_num in bins:
        print(f"Processing photo-z bin {bin_num}...")
        fn_bin = '/pscratch/sd/r/rhliu/projects/Weak_lensing/desi/photometric/' + f'desi_photo_LRG_bin{bin_num}.fits'
        cat_bin = Table.read(fn_bin)
        
        z_mean = redshifts[bin_num - 1]
        min_sep = physical_to_min_sep(physical_sep, z_mean, comoving=False, out_unit=u.arcmin)
        print(f"Photo-z bin {bin_num} mean z = {z_mean:.3f}, min_sep = {min_sep:.3f}")
        print('rp_to_arcmin:', rp_to_arcmin(1.0 * h, z_mean), '1 mpc to arcmin at z_mean')
        
        pruned = enforce_min_separation(cat_bin, ra_col="RA", dec_col="DEC", min_sep=min_sep)
        print(f"  Original count: {len(cat_bin)}, pruned count: {len(pruned)}")
        print('Pruned Fraction:', 1 - len(pruned)/len(cat_bin))
    pass