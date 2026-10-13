"""Truth tests for ``pyaraucaria.ffs``: compare FFS output with what was
injected into synthetic frames (``tests/ffs_synth.py``) or built by hand.
Production call sequences are guarded separately by ``test_ffs_regression``.
"""

import os
import tempfile
import unittest
import warnings

import numpy as np

from pyaraucaria.ffs import FFS
from tests.ffs_synth import FWHM_TO_SIGMA, gaussian_star, make_frame, scenario, star_field


def stamp(fwhm=4.0, ellipticity=0.0, theta=0.0, x=15.0, y=15.0, size=31, flux=1e5):
    img = np.zeros((size, size))
    sl, st = gaussian_star(img.shape, x, y, flux, fwhm, ellipticity, theta, radius=4 * size)
    img[sl] += st
    return img


def axial_diff(a, b):
    """Smallest difference between two orientations (period pi), radians."""
    return np.abs((a - b + np.pi / 2) % np.pi - np.pi / 2)


def match(found_xy, true_xy, tol=1.5):
    """Indices of true stars matched by a found star within ``tol`` px."""
    found_xy = np.asarray(found_xy, dtype=float).reshape(-1, 2)
    matched = set()
    extra = 0
    for fx, fy in found_xy:
        d = np.hypot(true_xy[:, 0] - fx, true_xy[:, 1] - fy)
        if d.size and d.min() <= tol:
            matched.add(int(d.argmin()))
        else:
            extra += 1
    return matched, extra


class QuietTestCase(unittest.TestCase):
    def setUp(self):
        self._w = warnings.catch_warnings()
        self._w.__enter__()
        warnings.simplefilter("ignore", RuntimeWarning)

    def tearDown(self):
        self._w.__exit__(None, None, None)


class TestStatic(QuietTestCase):
    """Per-stamp measurement helpers on noiseless stars."""

    def test_fwhm_round(self):
        fx, fy = FFS.fwhm(stamp(fwhm=4.0))
        self.assertAlmostEqual(fx, 4.0, delta=0.08)
        self.assertAlmostEqual(fy, 4.0, delta=0.08)

    def test_adaptive_moments_noiseless(self):
        for deg, ell, x, y in [(0, 0.0, 15, 15), (30, 0.3, 15.3, 14.6), (120, 0.4, 14.7, 15.2)]:
            with self.subTest(theta=deg, ellipticity=ell):
                img = stamp(fwhm=6.0, ellipticity=ell, theta=np.deg2rad(deg), x=x, y=y)
                f, e, t = FFS.adaptive_moments(img)
                self.assertAlmostEqual(e, ell, delta=1e-3)
                self.assertAlmostEqual(f, 6.0 * (1 - ell / 2), delta=0.01)
                if ell:
                    self.assertLess(axial_diff(t, np.deg2rad(deg)), np.deg2rad(0.1))

    def test_adaptive_moments_noise_unbiased(self):
        # Faint-ish star (peak SNR ~ 20) in a box with noise: mean e unbiased.
        rng = np.random.default_rng(0)
        clean = stamp(fwhm=5.0, ellipticity=0.3, theta=np.deg2rad(30), flux=2e4)
        es, ts = [], []
        for _ in range(50):
            f, e, t = FFS.adaptive_moments(clean + rng.normal(0, 22.0, clean.shape))
            es.append(e)
            ts.append(t)
        self.assertAlmostEqual(np.nanmean(es), 0.3, delta=0.03)
        self.assertLess(axial_diff(np.nanmedian(ts), np.deg2rad(30)), np.deg2rad(3))

    def test_adaptive_moments_failure_is_nan(self):
        self.assertTrue(np.isnan(FFS.adaptive_moments(np.zeros((20, 20)))[1]))
        self.assertTrue(np.isnan(FFS.adaptive_moments(-stamp())[1]))

    def test_concentration_index_gaussian(self):
        sigma = 4.0 * FWHM_TO_SIGMA
        analytic = (1 - np.exp(-4 / (2 * sigma ** 2))) / (1 - np.exp(-16 / (2 * sigma ** 2)))
        self.assertAlmostEqual(FFS.concentration_index(stamp(fwhm=4.0), 2, 4), analytic, delta=0.02)


class TestFrameStats(QuietTestCase):

    def test_background_and_noise(self):
        f = scenario("empty")      # sky 500 ADU, gain 1, read noise 5 e-
        ffs = FFS(f.image, gain=1.0, rn_noise=5.0)
        ffs.mk_stats()
        expected_sigma = np.sqrt(500 + 25)
        self.assertAlmostEqual(ffs.median, 500, delta=0.5)
        self.assertAlmostEqual(ffs.q_sigma, expected_sigma, delta=0.02 * expected_sigma)
        self.assertAlmostEqual(ffs.noise, expected_sigma, delta=0.02 * expected_sigma)
        self.assertAlmostEqual(ffs.rms, expected_sigma, delta=0.02 * expected_sigma)


class TestFindStars(QuietTestCase):

    def _find(self, f, **kw):
        ffs = FFS(f.image)
        ffs.mk_stats()
        ffs.find_stars(**kw)
        xy = np.column_stack([ffs.coo[:, 1], ffs.coo[:, 0]]) if len(ffs.coo) else np.empty((0, 2))
        return ffs, xy

    def test_clean_field_complete_and_pure(self):
        f = scenario("clean")
        _, xy = self._find(f, threshold=10, fwhm=4)
        truth = np.array([[s["x"], s["y"]] for s in f.stars])
        matched, extra = match(xy, truth, tol=1.0)
        self.assertEqual(len(matched), len(truth))
        self.assertEqual(extra, 0)

    def test_empty_field_no_detections(self):
        _, xy = self._find(scenario("empty"), threshold=5, fwhm=4)
        self.assertEqual(len(xy), 0)

    def test_sorted_by_peak(self):
        ffs, _ = self._find(scenario("clean"), threshold=10, fwhm=4)
        self.assertTrue(np.all(np.diff(ffs.adu) <= 0))

    def test_hot_pixels_detected_without_filters(self):
        # Documented legacy behaviour: plain find_stars keeps hot pixels.
        f = scenario("typical")
        _, xy = self._find(f, threshold=5, fwhm=4)
        hot = np.argwhere(f.truth["mask_hot"])[:, ::-1]
        matched, _ = match(xy, hot.astype(float), tol=1.0)
        self.assertEqual(len(matched), len(hot))

    def test_hot_pixels_rejected_by_filters(self):
        f = scenario("typical")
        _, xy = self._find(f, threshold=5, fwhm=4, max_concentration=0.5,
                           min_pixels_above_threshold=5)
        hot = np.argwhere(f.truth["mask_hot"])[:, ::-1]
        matched, _ = match(xy, hot.astype(float), tol=1.0)
        self.assertEqual(len(matched), 0)


class TestStarInfo(QuietTestCase):

    def test_frame_fwhm_round_stars(self):
        ffs = FFS(scenario("clean").image)
        ffs.calc_frame_fwhm(threshold=10, fwhm=4, box=10, N_stars=50, clip=4)
        self.assertAlmostEqual(ffs.frame_fwhm, 4.0, delta=0.2)
        self.assertLess(ffs.frame_ellipticity, 0.08)

    def test_rows_belong_to_detected_stars(self):
        f = scenario("clean")
        ffs = FFS(f.image)
        ffs.calc_frame_fwhm(threshold=10, fwhm=4, box=10, N_stars=50, clip=4)
        truth = np.array([[s["x"], s["y"]] for s in f.stars])
        matched, extra = match(np.column_stack([ffs.stars["x"], ffs.stars["y"]]), truth)
        self.assertEqual(extra, 0)
        self.assertEqual(len(matched), len(ffs.stars))

    def test_saturated_star_skipped_cleanly(self):
        # A skipped (saturated) star stays in the table with NaN measurements
        # (e.g. for the stars mask); the other rows are measured.
        stars = [dict(x=100.2, y=100.7, flux=5e6)]
        stars += [dict(x=60 + 60 * i + 0.3, y=300.4, flux=5e4) for i in range(6)]
        f = make_frame(shape=(400, 500), stars=stars, sky=500, saturation=45000, seed=3)
        ffs = FFS(f.image)
        ffs.saturation = 45000
        ffs.calc_frame_fwhm(threshold=10, fwhm=4, box=10, N_stars=50, clip=4)
        sat = np.asarray(ffs.stars["max_adu"]) >= 45000
        fwhm = np.asarray(ffs.stars["fwhm"])
        self.assertEqual(len(ffs.stars), 7)
        self.assertEqual(sat.sum(), 1)
        self.assertTrue(np.all(np.isnan(fwhm[sat])))
        self.assertTrue(np.all(np.isfinite(fwhm[~sat])))

    def test_frame_ellipticity(self):
        # Unweighted moments over the whole box with noise clipped at zero
        # gave 0.07 here; adaptive moments are needed.
        ffs = FFS(scenario("elliptical").image)
        ffs.calc_frame_fwhm(threshold=10, fwhm=5, box=15, N_stars=50, clip=4)
        self.assertAlmostEqual(ffs.frame_ellipticity, 0.3, delta=0.05)

    def test_frame_theta_spread(self):
        # theta has period pi; a plain std would explode on sign flips.
        ffs = FFS(scenario("elliptical").image)
        ffs.calc_frame_fwhm(threshold=10, fwhm=5, box=15, N_stars=50, clip=4)
        self.assertLess(ffs.frame_theta_spread, np.deg2rad(10))

    def test_frame_theta_value(self):
        ffs = FFS(scenario("elliptical").image)
        ffs.calc_frame_fwhm(threshold=10, fwhm=5, box=15, N_stars=50, clip=4)
        theta = np.asarray(ffs.theta)
        theta = theta[np.isfinite(theta)]
        med = 0.5 * np.angle(np.mean(np.exp(2j * theta)))   # axial mean
        self.assertLess(axial_diff(med, np.deg2rad(30)), np.deg2rad(5))


class TestSkyGradient(QuietTestCase):

    def test_surface_matches_model(self):
        f = scenario("gradient")
        ffs = FFS(f.image)
        ffs.sky_gradient()
        x, y = np.asarray(ffs.sky_surface_x), np.asarray(ffs.sky_surface_y)
        model = f.truth["sky_model"][y, x]
        np.testing.assert_allclose(ffs.sky_surface_bkg, model, atol=3.0)
        # amplitude over the whole frame, not only at segment centres
        full = f.truth["sky_model"]
        self.assertAlmostEqual(ffs.max_amplitude, full.max() - full.min(),
                               delta=0.02 * (full.max() - full.min()))

    def test_n_segments_honoured(self):
        ffs = FFS(scenario("gradient").image)
        ffs.sky_gradient(n_segments=5)
        self.assertEqual(len(ffs.sky_surface_x), 25)


class TestLines(QuietTestCase):

    def setUp(self):
        super().setUp()
        f = scenario("typical")   # one trail from (0, 60) to (511, 300)
        self.ffs = FFS(f.image)
        self.ffs.mk_stats()
        self.ffs.find_lines()
        d = np.array([511.0, 240.0])
        n = np.array([-d[1], d[0]]) / np.hypot(*d)
        self.theta_true = np.arctan2(n[1], n[0])
        self.rho_true = 60.0 * n[1]

    def _strongest(self):
        t, r = self.ffs.lines_theta[0], self.ffs.lines_rho[0]
        # (theta, rho) and (theta - pi, -rho) are the same line
        if abs(t - self.theta_true) > np.pi / 2:
            t, r = t + np.pi * np.sign(self.theta_true - t), -r
        return t, r

    def _offset_at_trail_centre(self):
        # Distance of the detected line from the middle of the true trail.
        # rho alone is not compared with the truth: it depends on theta, and
        # with a 1-deg theta grid a ~0.2 deg angle error shifts rho by ~1 px
        # (lever arm ~300 px) even when the line passes through the trail.
        t, r = self.ffs.lines_theta[0], self.ffs.lines_rho[0]
        xm, ym = 511.0 / 2, (60.0 + 300.0) / 2
        return abs(xm * np.cos(t) + ym * np.sin(t) - r)

    def test_strongest_line_orientation(self):
        t, _ = self._strongest()
        self.assertLess(abs(t - self.theta_true), np.deg2rad(1.5))
        self.assertLess(self._offset_at_trail_centre(), 1.5)

    def test_strongest_line_through_trail_centre(self):
        # The rho grid must agree with the voting rule (bin k <-> rho = k - rh0);
        # with linspace(-rh0, rh0, 2*rh0) the line was off by up to ~1 px.
        self.assertLess(self._offset_at_trail_centre(), 0.5)

    def test_single_trail_single_line(self):
        # The trail must be reported exactly once (the "butterfly" of partial
        # votes around it used to give extra lines). Other real lines in the
        # frame -- the hot column at x=260 -- are not counted here.
        xm, ym = 511.0 / 2, (60.0 + 300.0) / 2
        on_trail = [
            r for r, t in zip(self.ffs.lines_rho, self.ffs.lines_theta)
            if np.rad2deg(abs((t - self.theta_true + np.pi / 2) % np.pi - np.pi / 2)) < 3.0
            and abs(xm * np.cos(t) + ym * np.sin(t) - r) < 5.0
        ]
        self.assertEqual(len(on_trail), 1)


class TestExclude(QuietTestCase):
    """Sky map, NaN input and per-method exclusion masks."""

    def setUp(self):
        super().setUp()
        self.stars = star_field(n=60, seed=5)
        self.truth = np.array([[s["x"], s["y"]] for s in self.stars])
        self.image = make_frame(stars=self.stars, sky=1000, seed=6).image

    def _matched(self, ffs):
        return len(match(np.column_stack([ffs.coo[:, 1], ffs.coo[:, 0]]), self.truth, tol=2.0)[0])

    def test_sky_map_recovers_stars_on_gradient(self):
        xx = np.mgrid[0:self.image.shape[0], 0:self.image.shape[1]][1]
        ffs = FFS(self.image + 3.0 * xx)        # ~1500 ADU gradient across the frame
        ffs.mk_stats()
        ffs.find_stars(threshold=5, fwhm=4)     # sky map computed on demand
        self.assertIn("sky", ffs.maps)
        self.assertGreaterEqual(self._matched(ffs), 45)

    def test_nan_pixels_do_not_break_find_stars(self):
        img = self.image.copy()
        img[:, :20] = np.nan
        ffs = FFS(img)
        ffs.mk_stats()
        ffs.find_stars(threshold=5, fwhm=4)
        self.assertGreaterEqual(self._matched(ffs), 45)
        self.assertTrue(np.all(np.isfinite(ffs.adu)))

    def test_excluded_pixels_are_not_stars(self):
        img = self.image.copy()
        hot = [(100, 37), (300, 211), (450, 480)]
        mask = np.zeros(img.shape, bool)
        for y, x in hot:
            img[y, x] += 20000
            mask[y, x] = True
        ffs = FFS(img)
        ffs.masks["hot"] = mask
        ffs.mk_stats(exclude=FFS.EXCLUDE_ALL + ("hot",))
        ffs.find_stars(threshold=5, fwhm=4, exclude=FFS.EXCLUDE_DETECTION + ("hot",))
        self.assertFalse({tuple(c) for c in ffs.coo.tolist()} & set(hot))


class TestCalibration(QuietTestCase):
    """Bad columns / rows and hot / cold pixels on calibration frames (no stars)."""

    def setUp(self):
        super().setUp()
        self.rng = np.random.default_rng(1)

    @staticmethod
    def _found(ffs):
        cols = sorted(set(np.nonzero(ffs.masks["bad_columns"])[1].tolist()))
        rows = sorted(set(np.nonzero(ffs.masks["bad_rows"])[0].tolist()))
        return cols, rows

    def test_bias_columns_and_row(self):
        img = 1000 + self.rng.normal(0, 5, (512, 512))
        img[:, 100] += 10
        img[:, 300:302] -= 10       # 2 px wide, dark
        img[200, :] += 10
        ffs = FFS(np.round(img).astype(np.uint16))
        ffs.mk_columns_map()
        ffs.mk_bad_columns_mask(threshold=5)
        self.assertEqual(self._found(ffs), ([100, 300, 301], [200]))

    def test_partial_column_in_dark_with_cosmics(self):
        img = 1000 + self.rng.normal(0, 5, (512, 512))
        hits = self.rng.integers(0, 512, (500, 2))
        img[hits[:, 0], hits[:, 1]] += self.rng.uniform(100, 5000, 500)
        img[256:, 150] += 20        # defect from row 256 down
        ffs = FFS(img)
        ffs.mk_columns_map()
        ffs.mk_bad_columns_mask(threshold=5)
        self.assertEqual(self._found(ffs), ([150], []))
        self.assertTrue(ffs.masks["bad_columns"][256:, 150].all())
        self.assertFalse(ffs.masks["bad_columns"][:256, 150].any())

    def test_dark_hot_pixels(self):
        img = 1000 + self.rng.normal(0, 5, (256, 256))
        hot = self.rng.integers(0, 256, (50, 2))
        img[hot[:, 0], hot[:, 1]] += 200
        ffs = FFS(img)
        ffs.mk_pixels_map()
        ffs.mk_hot_pixels_mask(threshold=50)
        ffs.mk_cold_pixels_mask(threshold=50)
        truth = np.zeros(img.shape, bool)
        truth[hot[:, 0], hot[:, 1]] = True
        self.assertTrue(np.array_equal(ffs.masks["hot_pixels"], truth))
        self.assertFalse(ffs.masks["cold_pixels"].any())

    def test_flat_cold_pixels_with_vignetting(self):
        yy, xx = np.mgrid[0:256, 0:256]
        vign = 30000 * (1 - 0.3 * (np.hypot(xx - 128, yy - 128) / 128) ** 2)
        sens = np.ones((256, 256))
        cold = self.rng.integers(0, 256, (50, 2))
        sens[cold[:, 0], cold[:, 1]] = 0.8
        ffs = FFS(self.rng.poisson(vign * sens).astype(float))
        ffs.mk_pixels_map()
        ffs.mk_cold_pixels_mask(threshold=2000)
        truth = np.zeros((256, 256), bool)
        truth[cold[:, 0], cold[:, 1]] = True
        self.assertTrue(np.array_equal(ffs.masks["cold_pixels"], truth))


class TestLineMask(QuietTestCase):
    """line_to_mask: (rho, theta) and segment ends -> pixels."""

    shape = (100, 120)   # (ny, nx)

    def test_horizontal_full_line(self):
        # theta = pi/2: y = rho
        ys, xs = np.nonzero(FFS.line_to_mask(self.shape, rho=40, theta=np.pi / 2))
        self.assertTrue(np.all(ys == 40))
        self.assertEqual(xs.size, self.shape[1])

    def test_segment_is_limited_to_ends(self):
        m = FFS.line_to_mask(self.shape, rho=40, theta=np.pi / 2, half_width=0.5,
                             x0=10, y0=40, x1=50, y1=40)
        xs = np.nonzero(m)[1]
        self.assertEqual((xs.min(), xs.max()), (10, 50))


class TestRingTest(QuietTestCase):
    """Antipodal ring test used for Hough votes: lines vote, blobs do not."""

    def test_line_votes_blobs_do_not(self):
        yy, xx = np.mgrid[0:121, 0:121]
        line = np.abs((yy - 30) - 0.5 * (xx - 10)) <= 1.5       # ~27 deg, width ~3 px
        star = np.hypot(xx - 30, yy - 90) <= 3                   # smaller than the ring
        blob = np.hypot(xx - 95, yy - 20) <= 20                  # saturated star, larger than the ring
        votes = FFS._ring_test(line | star | blob, radius=8, omax=0.4)
        # the line votes along its whole length (edge pixels of a wide line may not)
        cols = np.nonzero((votes & line).any(axis=0))[0]
        self.assertTrue(set(range(20, 101)) <= set(cols.tolist()))
        self.assertFalse(votes[star | blob].any())


class TestSaveMef(QuietTestCase):

    def test_masks_and_tables_roundtrip(self):
        from astropy.io import fits
        ffs = FFS(scenario("typical").image, saturation=45000)
        ffs.calc_frame_fwhm(threshold=10, fwhm=4, box=10, N_stars=20)
        with tempfile.TemporaryDirectory() as d:
            path = os.path.join(d, "ffs.fits")
            ffs.save_mef(path)
            with fits.open(path) as hl:
                bits = hl["MASKS"].data
                for i, name in enumerate(ffs.masks):
                    self.assertEqual(hl["MASKS"].header[f"MASK{i}"], name)
                    np.testing.assert_array_equal((bits >> i) & 1 == 1, ffs.masks[name])
                self.assertEqual(len(hl["STARS"].data), len(ffs.stars))
                self.assertAlmostEqual(hl[0].header["FFFWHM"], ffs.frame_fwhm)
                self.assertEqual(hl["SKY"].data.shape, ffs.image.shape)


if __name__ == "__main__":
    unittest.main()
