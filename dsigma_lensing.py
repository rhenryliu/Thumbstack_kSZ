"""Weights of the 'DSigma_lensing' aperture filter: the lensing-bin average of DeltaSigma.

A weak-lensing DeltaSigma measurement in a radial bin [a, b] reports the area-weighted mean of
DeltaSigma(R) = mean Sigma(<R) - Sigma(R) over the bin. That mean is a compensated aperture filter
on Sigma,

    <DS>_[a,b] = int d^2r Sigma(r) w(r),
    w(r) = [2 ln(b/a) (r < a);  2 ln(b/r) - 1 (a <= r < b);  0 (r >= b)] / (pi (b^2 - a^2)),

so applying w to a kSZ stamp measures the gas exactly as lensing measures the total mass, in the
same bins. The bins are the edges whose geometric means are the aperture radii, as the DESI x HSC
lensing bins are built (DESIxHSC-Lensing docs/ksz_lensing_matched_filter.md).

bins_from_geomean_monotonic, _signed_quadrant_area, pixel_disk_overlap and
delta_sigma_lensing_weights are copied from SimulationStacker (src/utils.py and src/filters.py,
commit da39be7), which uses the same weights as filterType 'DSigma_lensing'. The repositories are
kept independent on purpose, so this is a copy rather than an import; keep the two in step.
gaussian_lensing_bin_profile is new here: a vectorised form of the exact Gaussian bin average in
SimulationStacker's tests (tests/test_delta_sigma_lensing.py), times pi R0^2.
"""
from typing import Tuple

import numpy as np
from scipy.special import exp1


def bins_from_geomean_monotonic(r_desired, bins0=None, pick='mid'):
    """Recover strictly increasing bin edges from desired geometric-mean centres.

    Given a positive, strictly increasing array r_desired of length N, where each element is the
    geometric mean of consecutive edges (r_desired[i] = sqrt(bins[i] * bins[i+1])), return the
    N + 1 edges. The leftmost edge bins0 must lie in an open interval (lower, upper) for the
    edges to increase; it is validated if given, otherwise chosen by pick.

    Args:
      r_desired (array_like): Length-N positive array of target geometric means.
      bins0 (float, optional): Leftmost edge; must satisfy lower < bins0 < upper.
      pick (str): How to choose bins0 if not given: 'mid'/'mean' (arithmetic midpoint of
        (lower, upper)), 'geom'/'geometric' (geometric mean), 'lower' or 'upper'.

    Returns:
      Tuple[np.ndarray, Tuple[float, float]]: the (N+1,) edges and the interval (lower, upper).

    Raises:
      ValueError: If r_desired has a non-positive element, no increasing solution exists,
        bins0 is outside (lower, upper), or pick is not recognised.
      RuntimeError: If the constructed edges are not strictly increasing.
    """
    r_desired = np.asarray(r_desired, dtype=float)
    if np.any(r_desired <= 0):
        raise ValueError("r_desired must be positive")
    n = len(r_desired)

    # allowed interval for bins0, propagating k_i with alternating exponent sign
    lower = 0.0
    upper = np.inf
    k = 1.0
    s = 1  # +1: even index b_i = k*b0; -1: odd index b_i = k/b0
    for i in range(n):
        if s == 1:
            upper = min(upper, r_desired[i] / k)      # k * b0 < r_desired[i]
        else:
            lower = max(lower, k / r_desired[i])      # k / b0 < r_desired[i]
        k = r_desired[i]**2 / k
        s = -s

    if not (lower < upper and upper > 0):
        raise ValueError("No monotonic increasing solution for given r_desired")

    if bins0 is None:
        if pick in ('mid', 'mean'):
            bins0 = 0.5 * (lower + upper)
        elif pick in ('geom', 'geometric'):
            bins0 = np.sqrt(lower * upper)
        elif pick == 'lower':
            bins0 = lower + 1e-12
        elif pick == 'upper':
            bins0 = upper - 1e-12
        else:
            raise ValueError("unknown pick option")
    else:
        bins0 = float(bins0)
        if not (lower < bins0 < upper):
            raise ValueError(f"bins0={bins0} outside allowed interval ({lower}, {upper})")

    bins = np.empty(n + 1, dtype=float)
    bins[0] = bins0
    for i in range(n):
        bins[i + 1] = r_desired[i]**2 / bins[i]

    if not np.all(bins[1:] > bins[:-1]):
        raise RuntimeError("constructed bins are not strictly increasing (unexpected)")

    return bins, (lower, upper)


def _signed_quadrant_area(x, y, R):
    """Signed area of the disk u^2 + v^2 < R^2 between the origin and the corner (x, y).

    For x, y >= 0 this is the area of the disk inside [0, x] x [0, y]. Flipping the sign of
    x or y flips the sign of the result, so the area of the disk inside any axis-aligned
    rectangle [x0, x1] x [y0, y1] is Q(x1, y1) - Q(x0, y1) - Q(x1, y0) + Q(x0, y0).

    Args:
        x (np.ndarray): Corner abscissae (broadcastable with ``y`` and ``R``).
        y (np.ndarray): Corner ordinates.
        R (np.ndarray): Disk radii, > 0.

    Returns:
        np.ndarray: Signed areas, in the squared units of the inputs.
    """
    sign = np.sign(x) * np.sign(y)
    x = np.minimum(np.abs(x), R)
    y = np.minimum(np.abs(y), R)
    u_star = np.sqrt(np.maximum(R**2 - y**2, 0.0))   # where the circle drops to height y

    def chord_integral(u):
        """Integral of sqrt(R^2 - t^2) dt from 0 to u."""
        return 0.5 * (u * np.sqrt(np.maximum(R**2 - u**2, 0.0))
                      + R**2 * np.arcsin(np.clip(u / R, -1.0, 1.0)))

    area = np.where(x <= u_star, x * y, y * u_star + chord_integral(x) - chord_integral(u_star))
    return sign * area


def pixel_disk_overlap(x0, x1, y0, y1, R):
    """Exact area of the disk of radius R centred on the origin inside each rectangle.

    Args:
        x0 (np.ndarray): Lower row-axis edges of the rectangles.
        x1 (np.ndarray): Upper row-axis edges.
        y0 (np.ndarray): Lower column-axis edges.
        y1 (np.ndarray): Upper column-axis edges.
        R (np.ndarray): Disk radii, > 0, broadcastable with the edges.

    Returns:
        np.ndarray: Overlap areas, in the squared units of the inputs.
    """
    Q = _signed_quadrant_area
    return Q(x1, y1, R) - Q(x0, y1, R) - Q(x1, y0, R) + Q(x0, y0, R)


def delta_sigma_lensing_weights(half: int, pixel_size: float, edges: np.ndarray,
                                offset: Tuple[float, float] = (0.0, 0.0),
                                n_quad: int = 256) -> np.ndarray:
    """Per-pixel weights that turn a cutout into the lensing-bin average of DeltaSigma.

    For a bin [a, b], with M(<R) the projected mass (or flux) within R of the centre,

        <DS>_[a,b] = [int_a^b 2 M(<R) / R dR - (M(<b) - M(<a))] / (pi (b^2 - a^2)).

    Treating the content of each pixel as uniform over the pixel, M(<R) = sum_p m_p f_p(R), where
    f_p(R) is the exact fraction of pixel p inside radius R (pixel_disk_overlap), so
    <DS> = sum_p m_p W_p with

        W_p = [int_a^b 2 f_p(R) / R dR - f_p(b) + f_p(a)] / (pi (b^2 - a^2)).

    The R integral uses Gauss-Legendre quadrature. Since sum_p f_p(R) x (pixel area) = pi R^2
    exactly once the cutout contains the disk of radius b, the weights cancel a uniform sheet to
    rounding error for any number of nodes.

    Args:
        half: Cutout half-width in pixels; the cutout is (2*half+1, 2*half+1), centred on its
            central pixel.
        pixel_size: Pixel size, in the units of ``edges``.
        edges: Bin edges, shape (N+1,), positive and strictly increasing.
        offset: (row, column) position of the centre relative to the centre of the central
            pixel, in the units of ``edges``. ThumbStack's stamps are centred on the object, so
            (0, 0) there.
        n_quad: Gauss-Legendre nodes per bin.

    Returns:
        Weights of shape (N, 2*half+1, 2*half+1), per unit area of ``edges``. For a cutout ``m``
        of content per pixel, ``np.tensordot(W, m, axes=([1, 2], [0, 1]))`` is <DS> in each bin;
        for a cutout of surface density (map units), multiply by the pixel area.

    Raises:
        ValueError: If the edges are not positive and increasing, or if the disk of radius
            ``edges[-1]`` around the centre is not inside the cutout.
    """
    edges = np.asarray(edges, dtype=float)
    if edges[0] <= 0 or np.any(np.diff(edges) <= 0):
        raise ValueError(f'edges must be positive and strictly increasing, got {edges}')
    off = np.asarray(offset, dtype=float)
    reach = (half + 0.5) * pixel_size - np.max(np.abs(off))
    if reach < edges[-1]:
        raise ValueError(f'cutout reaches {reach:.4g} from the centre, inside the outer bin edge '
                         f'{edges[-1]:.4g}; enlarge the cutout')

    n = 2 * half + 1
    centres = np.arange(-half, half + 1) * pixel_size
    # Pixel rectangles relative to the centre, flattened row-major.
    x0 = np.repeat(centres - 0.5 * pixel_size - off[0], n)
    y0 = np.tile(centres - 0.5 * pixel_size - off[1], n)
    x1, y1 = x0 + pixel_size, y0 + pixel_size
    nearest = lambda lo, hi: np.where((lo <= 0) & (hi >= 0), 0.0, np.minimum(np.abs(lo), np.abs(hi)))
    r_near = np.hypot(nearest(x0, x1), nearest(y0, y1))
    r_far = np.hypot(np.maximum(np.abs(x0), np.abs(x1)), np.maximum(np.abs(y0), np.abs(y1)))
    area = pixel_size**2

    nodes, node_weights = np.polynomial.legendre.leggauss(n_quad)
    W = np.zeros((len(edges) - 1, n * n))
    for i, (a, b) in enumerate(zip(edges[:-1], edges[1:])):
        R = 0.5 * (b - a) * nodes + 0.5 * (b + a)
        wq = 0.5 * (b - a) * node_weights * 2.0 / R          # quadrature weights of int 2 (.) / R dR
        inside = r_far <= a                                  # f_p = 1 on [a, b]
        W[i, inside] = wq.sum()                              # (int 2/R dR) - 1 + 1
        p = np.nonzero((r_near < b) & ~inside)[0]            # f_p changes somewhere on [a, b]
        if p.size:
            frac = lambda r: pixel_disk_overlap(x0[p], x1[p], y0[p], y1[p], r) / area
            W[i, p] = wq @ frac(R[:, None]) - frac(b) + frac(a)
        W[i] /= np.pi * (b**2 - a**2)
    return W.reshape(len(edges) - 1, n, n)


def gaussian_lensing_bin_profile(edges, r_centres, sigma):
    """pi R0^2 <DS>_[a,b] for a unit-flux circular Gaussian of width sigma (exact).

    With x = R^2 / (2 sigma^2), int_a^b DS(R) 2R dR = [Ein(x) + e^-x] / pi between the edges,
    where Ein(x) = E1(x) + ln x + gamma. Multiplying the bin average by pi R0^2 makes the result
    dimensionless; it tends to 2 R0^2 ln(b/a) / (b^2 - a^2) (about 1) for a point source.

    Args:
        edges (np.ndarray): Bin edges, shape (N+1,), same units as sigma.
        r_centres (np.ndarray): Aperture radii R0, shape (N,).
        sigma (float): Gaussian width.

    Returns:
        np.ndarray: Shape (N,).
    """
    edges = np.asarray(edges, dtype=float)

    def F(r):
        x = r**2 / (2.0 * sigma**2)
        series = x - x**2 / 4 + x**3 / 18 - x**4 / 96          # Ein(x) for small x
        ein = np.where(x < 1e-3, series, exp1(np.maximum(x, 1e-300)) + np.log(np.maximum(x, 1e-300))
                       + np.euler_gamma)
        return ein + np.exp(-x)

    a, b = edges[:-1], edges[1:]
    return np.asarray(r_centres)**2 * (F(b) - F(a)) / (b**2 - a**2)
