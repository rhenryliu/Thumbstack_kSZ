"""Tests of the 'DSigma_lensing' filter (dsigma_lensing.py and its branches in thumbstack.py).

Known answers use a circular Gaussian of width s, whose lensing-bin mean of DeltaSigma over
[a, b] is exact: with x = R^2 / 2 s^2 and F(x) = Ein(x) + e^-x,
<DS>_[a,b] = [F(x_b) - F(x_a)] / (pi (b^2 - a^2)) for unit flux.
Run in the rhliu_tSZ environment from the repository root: python -m pytest tests/
"""
import os
import sys
import types

import numpy as np
import pytest

sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..'))

import dsigma_lensing as dl  # noqa: E402

R0S = np.linspace(1.0, 6.0, 9)
EDGES = dl.bins_from_geomean_monotonic(R0S, pick='geom')[0]
ARCMIN = np.pi / 180. / 60.


def exact_bin_mean(a, b, s):
    """Lensing-bin mean of DeltaSigma for a unit-flux Gaussian (pi R0^2 removed)."""
    r0 = np.sqrt(a * b)
    return dl.gaussian_lensing_bin_profile(np.array([a, b]), np.array([r0]), s)[0] / (np.pi * r0**2)


def gaussian(r, s):
    return np.exp(-r**2 / (2 * s**2)) / (2 * np.pi * s**2)


def test_edges_are_the_lensing_bins():
    """Geometric-mean edges of the 1-6' grid, as the DESI x HSC lensing bins."""
    np.testing.assert_allclose(np.sqrt(EDGES[:-1] * EDGES[1:]), R0S, rtol=1e-12)
    np.testing.assert_allclose(EDGES[[0, -1]], [0.74056969, 6.33924653], rtol=1e-7)


def test_uniform_sheet_cancels_and_point_mass_is_exact():
    half, pix = 34, 0.25
    W = dl.delta_sigma_lensing_weights(half, pix, EDGES)
    scale = np.abs(W).sum(axis=(1, 2)) * pix**2
    assert np.all(np.abs(W.sum(axis=(1, 2)) * pix**2) <= 1e-12 * scale)
    m = np.zeros((2 * half + 1,) * 2)
    m[half, half] = 1.0
    a, b = EDGES[:-1], EDGES[1:]
    np.testing.assert_allclose(np.tensordot(W, m, axes=([1, 2], [0, 1])),
                               2 * np.log(b / a) / (np.pi * (b**2 - a**2)), rtol=1e-12)


@pytest.mark.parametrize('fwhm', [1.7, 2.3, 4.0])
def test_gaussian_on_ideal_stamp_converges(fwhm):
    """Gaussian sampled at the stamp's pixel centres (as ThumbStack's stamps are): within 3% on
    the production 0.25' stamp, and closer on a 0.125' stamp."""
    s = fwhm / np.sqrt(8 * np.log(2))
    exact = np.array([exact_bin_mean(a, b, s) for a, b in zip(EDGES[:-1], EDGES[1:])])
    err = []
    for pix in (0.25, 0.125):
        half = int(np.ceil(8.625 / pix))
        c = np.arange(-half, half + 1) * pix
        stamp = gaussian(np.hypot(c[:, None], c[None, :]), s)
        W = dl.delta_sigma_lensing_weights(half, pix, EDGES)
        got = np.tensordot(W, stamp * pix**2, axes=([1, 2], [0, 1]))
        err.append(np.max(np.abs(got / exact - 1)))
    assert err[0] < 0.03
    assert err[1] < err[0] / 2


def test_theory_point_source_limit():
    got = dl.gaussian_lensing_bin_profile(EDGES, R0S, 1e-3)
    a, b = EDGES[:-1], EDGES[1:]
    np.testing.assert_allclose(got, 2 * R0S**2 * np.log(b / a) / (b**2 - a**2), rtol=1e-10)


@pytest.fixture(scope='module')
def thumbstack():
    """thumbstack needs the rhliu_tSZ environment; skip only the tests that use it."""
    return pytest.importorskip('thumbstack')


def test_aperture_grid_matches_names_exactly(thumbstack):
    """After the fix at __init__'s call, loadAPRadii gets the parsed array, so names are matched
    exactly: 'DSigma_lensing' takes the 1-6' grid, today's names are unchanged, and a name that
    only contains 'DSigma' as a substring matches nothing."""
    load = thumbstack.ThumbStack.loadAPRadii
    for names in (['DSigma'], ['DSigma_lensing'], ['diskring'], ['upsilon']):
        st = types.SimpleNamespace()
        load(st, np.array(names))
        np.testing.assert_allclose(st.RApArcmin, R0S)
    st = types.SimpleNamespace()
    load(st, np.array(['meanring']))
    np.testing.assert_allclose(st.RApArcminBins, np.linspace(0, 10, 10))
    st = types.SimpleNamespace()
    load(st, np.array(['xDSigmax']))
    assert not hasattr(st, 'RApArcmin')
    st = types.SimpleNamespace()
    load(st, np.array(['DSigma_lensing']), rApMinArcmin=1.0, rApMaxArcmin=9.75, nRAp=15)
    np.testing.assert_allclose(st.RApArcmin, np.linspace(1.0, 9.75, 15))


def test_analyze_object_end_to_end(thumbstack):
    """'DSigma_lensing' through ThumbStack's own extractStamp and analyzeObject equals the
    weights applied directly to the same stamp, and the stamp size holds the outer edge."""
    from pixell import enmap
    TS = thumbstack.ThumbStack
    st = types.SimpleNamespace(resCutoutArcmin=0.25, projCutout='cea', cmbHit=None,
                               filterTypes=np.array(['DSigma_lensing']), overlapFlag=np.array([1]))
    TS.loadAPRadii(st, st.filterTypes)
    for name in ('cutoutGeometry', 'extractStamp', 'aperturePhotometryFilter', 'loadDSigmaLensingWeights'):
        setattr(st, name, types.MethodType(getattr(TS, name), st))
    st.loadDSigmaLensingWeights()
    assert st.DSigmaLensingWeights.shape == (9, 69, 69)

    # 0.5' CAR map with a 2.3'-FWHM Gaussian at a position off the pixel grid
    shape, wcs = enmap.geometry(pos=np.array([[-0.5, -0.5], [0.5, 0.5]]) * np.pi / 180, res=0.5 * ARCMIN, proj='car')
    ra0, dec0 = 0.0123, -0.0071                                     # deg
    pos = enmap.posmap(shape, wcs)
    r = np.hypot((pos[1] - ra0 * np.pi / 180) * np.cos(dec0 * np.pi / 180), pos[0] - dec0 * np.pi / 180) / ARCMIN
    s = 2.3 / np.sqrt(8 * np.log(2))
    st.cmbMap = enmap.ndmap(gaussian(r, s), wcs)
    st.cmbMask = enmap.ones(shape, wcs)
    st.Catalog = types.SimpleNamespace(RA=np.array([ra0]), DEC=np.array([dec0]), Z=np.array([0.5]))

    filt = TS.analyzeObject(st, 0)[0]['DSigma_lensing']
    opos, stamp, _, _ = st.extractStamp(ra0, dec0)
    area = float(opos[1].area() / opos[1].size) / ARCMIN**2
    # the weights assume the object sits on the centre of the central pixel, of area res^2
    np.testing.assert_allclose(opos[:, 34, 34], 0.0, atol=1e-12)
    np.testing.assert_allclose(area, st.resCutoutArcmin**2, rtol=1e-8)
    direct = np.tensordot(st.DSigmaLensingWeights, np.asarray(stamp) * area, axes=([1, 2], [0, 1]))
    np.testing.assert_allclose(filt, direct, rtol=1e-10)
    exact = np.array([exact_bin_mean(a, b, s) for a, b in zip(EDGES[:-1], EDGES[1:])])
    # the map's 0.5' pixels and the interpolation onto the stamp smooth the profile; for this
    # Gaussian the pipeline reads -7.3%, -5.2%, -3.0%, -1.3% at 1-2.9' and < 0.4% beyond
    assert np.all(np.abs(filt / exact - 1) < np.array([0.10, 0.08, 0.05, 0.03] + [0.02] * 5))
    np.testing.assert_allclose(st.RApEdgesArcmin, EDGES, rtol=1e-12)
    # the theory profile is the exact bin mean per unit integrated flux, per sr
    theory = TS.ftheoryGaussianProfile(st, s, filterType='DSigma_lensing')
    np.testing.assert_allclose(theory * ARCMIN**2, exact, rtol=1e-12)
