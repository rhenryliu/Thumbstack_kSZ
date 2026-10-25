"""Tests of the completeness-weighted kSZ estimator 'ksz_compweight' (thumbstack.py) and of the
weight reader desi/spec_Y3/weight_utils.py.

The estimator is sigma_w / rV * sum w v' T / sum w v'^2, with v' = v - <v>_w over the stacked
objects (after any resampling or shuffling) and sigma_w^2 = sum w v'^2 / sum w; the reference
below is written independently of thumbstack.py.
Run in the rhliu_tSZ environment from the repository root: python -m pytest tests/
"""
import os
import sys
import types

import numpy as np
import pandas as pd
import pytest

ROOT = os.path.join(os.path.dirname(__file__), '..')
sys.path.insert(0, ROOT)
sys.path.insert(0, os.path.join(ROOT, 'desi', 'spec_Y3'))

import weight_utils  # noqa: E402

FT = 'DSigma_lensing'
NRAP = 9


@pytest.fixture(scope='module')
def thumbstack():
    """thumbstack needs the rhliu_tSZ environment; skip only the tests that use it."""
    return pytest.importorskip('thumbstack')


def reference(t, vR, w):
    """Weighted estimator from its definition (v = -vR / c, as in ThumbStack)."""
    v = -vR / 3.e5
    vp = v - np.sum(w * v) / np.sum(w)
    sw = np.sqrt(np.sum(w * vp**2) / np.sum(w))
    return sw * np.sum((w * vp)[:, None] * t, axis=0) / np.sum(w * vp**2)


def make_catalogue(n=500, seed=1, weights='random'):
    """Temperatures with a kSZ-like term, a mean (tSZ-like) offset and noise; velocities with a
    non-zero mean, as for the stacked objects after zeroV and the mask."""
    rng = np.random.default_rng(seed)
    vR = rng.normal(20., 300., n)
    t = 0.3 + 0.05 * (-vR / 300.)[:, None] * np.linspace(1., 2., NRAP) + rng.normal(0., 1., (n, NRAP))
    if weights == 'random':
        w = rng.uniform(0.5, 3., n)
    elif weights == 'integer':
        w = rng.integers(1, 5, n).astype(float)
    else:
        w = np.ones(n)
    mask = rng.uniform(size=n) > 0.2
    return t, vR, w, mask


def stub(t, vR, w=None):
    """The attributes computeStackedProfile reads, with ts=None (ts is self)."""
    n = len(vR)
    cat = types.SimpleNamespace(vR=vR.copy(), Mvir=np.zeros(n), rV=1., integratedY=np.zeros(n),
                                integratedKSZ=np.zeros(n), RA=np.zeros(n), Z=np.full(n, 0.5))
    if w is not None:
        cat.wComp = w
    return types.SimpleNamespace(Catalog=cat, filtMap={FT: t}, filtVarTrue={FT: np.ones_like(t)},
                                 filtHitNoiseStdDev={FT: np.ones_like(t)}, mMin=0., mMax=np.inf,
                                 cmbUnitLatex='', use_median=False, nRAp=NRAP)


def stack(thumbstack, st, est, mask, **kwargs):
    return thumbstack.ThumbStack.computeStackedProfile(st, FT, est, mask=mask, **kwargs)[0]


def test_matches_reference(thumbstack):
    t, vR, w, mask = make_catalogue()
    got = stack(thumbstack, stub(t, vR, w), 'ksz_compweight', mask)
    np.testing.assert_allclose(got, reference(t[mask], vR[mask], w[mask]), rtol=1e-12)


def test_unit_weights_equal_uniform_recentred_after_mask(thumbstack):
    t, vR, _, mask = make_catalogue(weights='unit')
    got = stack(thumbstack, stub(t, vR, np.ones(len(vR))), 'ksz_compweight', mask)
    want = stack(thumbstack, stub(t, vR - vR[mask].mean()), 'ksz_uniformweight', mask)
    np.testing.assert_allclose(got, want, rtol=1e-12)


def test_integer_weights_equal_repeated_objects(thumbstack):
    t, vR, w, mask = make_catalogue(weights='integer')
    got = stack(thumbstack, stub(t, vR, w), 'ksz_compweight', mask)
    rep = np.repeat(np.arange(len(vR))[mask], w[mask].astype(int))
    want = stack(thumbstack, stub(t[rep], vR[rep], np.ones(len(rep))), 'ksz_compweight', np.ones(len(rep), bool))
    np.testing.assert_allclose(got, want, rtol=1e-12)


def test_constant_temperature_gives_zero(thumbstack):
    """sum w v' = 0: a temperature common to all objects (tSZ, dust, a monopole) cancels."""
    _, vR, w, mask = make_catalogue()
    t = np.full((len(vR), NRAP), 0.3)
    got = stack(thumbstack, stub(t, vR, w), 'ksz_compweight', mask)
    assert np.all(np.abs(got) < 1e-13)
    # the unweighted estimator leaks it, since the stacked velocities do not average to zero
    assert np.all(np.abs(stack(thumbstack, stub(t, vR), 'ksz_uniformweight', mask)) > 1e-3)


def test_velocity_offset_changes_nothing(thumbstack):
    """So zeroV, applied before the mask, does not affect ksz_compweight."""
    t, vR, w, mask = make_catalogue()
    a = stack(thumbstack, stub(t, vR, w), 'ksz_compweight', mask)
    b = stack(thumbstack, stub(t, vR - 123.4, w), 'ksz_compweight', mask)
    np.testing.assert_allclose(a, b, rtol=1e-10)


def test_bootstrap_resamples_weights_with_objects(thumbstack):
    t, vR, w, mask = make_catalogue()
    got = stack(thumbstack, stub(t, vR, w), 'ksz_compweight', mask, iBootstrap=7)
    np.random.seed(7)
    J = np.random.choice(np.arange(mask.sum()), size=mask.sum(), replace=True)
    np.testing.assert_allclose(got, reference(t[mask][J], vR[mask][J], w[mask][J]), rtol=1e-12)


def test_shuffle_keeps_weights_on_galaxies(thumbstack):
    t, vR, w, mask = make_catalogue()
    got = stack(thumbstack, stub(t, vR, w), 'ksz_compweight', mask, iVShuffle=3)
    np.random.seed(3)
    J = np.random.permutation(np.arange(mask.sum()))
    np.testing.assert_allclose(got, reference(t[mask], vR[mask][J], w[mask]), rtol=1e-12)


def test_catalogue_mask_drops_objects_without_weight(thumbstack):
    t, vR, w, _ = make_catalogue()
    st = stub(t, vR, w)
    st.overlapFlag = np.ones(len(vR))
    st.filtMask = {FT: np.zeros_like(t)}
    mask_all = thumbstack.ThumbStack.catalogMask(st, filterType=FT)
    del st.Catalog.wComp
    np.testing.assert_array_equal(mask_all, thumbstack.ThumbStack.catalogMask(st, filterType=FT))
    # without outlier rejection, whose threshold depends on the objects kept
    mask_all = thumbstack.ThumbStack.catalogMask(st, filterType=FT, outlierReject=False)
    st.Catalog.wComp = w.copy()
    st.Catalog.wComp[:10] = np.nan
    mask_nan = thumbstack.ThumbStack.catalogMask(st, filterType=FT, outlierReject=False)
    assert not mask_nan[:10].any()
    np.testing.assert_array_equal(mask_nan[10:], mask_all[10:])


# --- weight_utils ---

def write_table(path, n=50, seed=2):
    rng = np.random.default_rng(seed)
    df = pd.DataFrame({'TARGETID': np.arange(1000, 1000 + n), 'RA': rng.uniform(0, 360, n),
                       'DEC': rng.uniform(-10, 30, n), 'Z': rng.uniform(0.1, 0.4, n),
                       'WEIGHT_COMP': rng.integers(1, 4, n).astype(float),
                       'WEIGHT_ZFAIL': rng.uniform(1., 1.05, n)})
    df.to_csv(path, index=False)
    cat = types.SimpleNamespace(RA=df['RA'].values, DEC=df['DEC'].values, Z=df['Z'].values, nObj=n)
    return df, cat


def test_load_weights_product_of_columns(tmp_path):
    df, cat = write_table(tmp_path / 'cat.csv')
    w = weight_utils.load_weights({'file': str(tmp_path / 'cat.csv')}, cat, 'DESIY3_BGS', 'BGS', '', 'full', 'x')
    np.testing.assert_allclose(w, df['WEIGHT_COMP'].values * df['WEIGHT_ZFAIL'].values, rtol=1e-15)  # CSV parsing: 1 ulp
    w1 = weight_utils.load_weights({'file': str(tmp_path / 'cat.csv'), 'columns': ['WEIGHT_COMP']},
                                   cat, 'DESIY3_BGS', 'BGS', '', 'full', 'x')
    np.testing.assert_array_equal(w1, df['WEIGHT_COMP'].values)


def test_load_weights_default_path_is_the_intermediate_csv(tmp_path):
    df, cat = write_table(tmp_path / 'full_zbin2_x.csv')
    w = weight_utils.load_weights({}, cat, 'DESIY3_LRG_z2', 'LRG', str(tmp_path), 'full', 'x')
    np.testing.assert_allclose(w, df['WEIGHT_COMP'].values * df['WEIGHT_ZFAIL'].values, rtol=1e-15)  # CSV parsing: 1 ulp


def test_load_weights_refuses_misaligned_tables(tmp_path):
    _, cat = write_table(tmp_path / 'cat.csv')
    cfg = {'file': str(tmp_path / 'cat.csv')}
    shifted = types.SimpleNamespace(RA=np.roll(cat.RA, 1), DEC=cat.DEC, Z=cat.Z, nObj=cat.nObj)
    with pytest.raises(ValueError, match='RA'):
        weight_utils.load_weights(cfg, shifted, 'DESIY3_BGS', 'BGS', '', 'full', 'x')
    short = types.SimpleNamespace(RA=cat.RA[:-1], DEC=cat.DEC[:-1], Z=cat.Z[:-1], nObj=cat.nObj - 1)
    with pytest.raises(ValueError, match='rows'):
        weight_utils.load_weights(cfg, short, 'DESIY3_BGS', 'BGS', '', 'full', 'x')


def test_load_weights_refuses_nonpositive_weights(tmp_path):
    df, cat = write_table(tmp_path / 'cat.csv')
    df.loc[3, 'WEIGHT_COMP'] = 0.
    df.to_csv(tmp_path / 'cat.csv', index=False)
    with pytest.raises(ValueError, match='finite and positive'):
        weight_utils.load_weights({'file': str(tmp_path / 'cat.csv')}, cat, 'DESIY3_BGS', 'BGS', '', 'full', 'x')


def test_join_by_targetid_with_missing_objects(tmp_path):
    df, cat = write_table(tmp_path / 'cat.csv')
    other = df.sample(frac=1., random_state=0).iloc[5:]  # shuffled, 5 objects missing
    other = other.assign(WEIGHT_COMP=other['WEIGHT_COMP'] + 1.)
    other.iloc[:20].to_csv(tmp_path / 'ngc.csv', index=False)
    other.iloc[20:].to_csv(tmp_path / 'sgc.csv', index=False)
    cfg = {'file': str(tmp_path / 'cat.csv'), 'join': [str(tmp_path / 'ngc.csv'), str(tmp_path / 'sgc.csv')]}
    w = weight_utils.load_weights(cfg, cat, 'DESIY3_BGS', 'BGS', '', 'full', 'x')
    missing = ~df['TARGETID'].isin(other['TARGETID']).values
    assert missing.sum() == 5 and np.all(np.isnan(w[missing]))
    want = (df['WEIGHT_COMP'].values + 1.) * df['WEIGHT_ZFAIL'].values
    np.testing.assert_allclose(w[~missing], want[~missing], rtol=1e-15)
    other.iloc[:2].to_csv(tmp_path / 'dup.csv', index=False)
    with pytest.raises(ValueError, match='duplicate'):
        weight_utils.load_weights({**cfg, 'join': cfg['join'] + [str(tmp_path / 'dup.csv')]},
                                  cat, 'DESIY3_BGS', 'BGS', '', 'full', 'x')
    other.assign(TARGETID=other['TARGETID'] + 10**6).to_csv(tmp_path / 'none.csv', index=False)
    with pytest.raises(ValueError, match='no object has a weight'):
        weight_utils.load_weights({**cfg, 'join': [str(tmp_path / 'none.csv')]}, cat, 'DESIY3_BGS', 'BGS', '', 'full', 'x')


def test_join_keeps_large_targetids_distinct(tmp_path):
    """DESI TARGETIDs exceed 2^53: a float64 comparison would match neighbouring IDs."""
    df, cat = write_table(tmp_path / 'cat.csv')
    big = 39627000000000001 + np.arange(len(df))
    df.assign(TARGETID=big).to_csv(tmp_path / 'cat.csv', index=False)
    df.assign(TARGETID=big + 1).iloc[:-1].to_csv(tmp_path / 'shifted.csv', index=False)
    cfg = {'file': str(tmp_path / 'cat.csv'), 'join': [str(tmp_path / 'shifted.csv')]}
    w = weight_utils.load_weights(cfg, cat, 'DESIY3_BGS', 'BGS', '', 'full', 'x')
    assert np.isnan(w[0]) and np.all(np.isfinite(w[1:]))
    np.testing.assert_allclose(w[1:], (df['WEIGHT_COMP'] * df['WEIGHT_ZFAIL']).values[:-1], rtol=1e-15)
