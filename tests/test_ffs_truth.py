"""Truth tests for ``pyaraucaria.ffs``: compare FFS output with what was
injected into synthetic frames (``tests/ffs_synth.py``).

Known FFS defects are kept as ``@unittest.expectedFailure`` tests labelled
``KNOWN BUG``. They document the defect without changing production
behaviour; when a fix lands the test reports "unexpected success" and the
decorator should be removed (and the regression snapshot regenerated).
"""

import unittest
import warnings

import numpy as np

from pyaraucaria.ffs import FFS
from tests.ffs_synth import FWHM_TO_SIGMA, gaussian_star, make_frame, scenario


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


class TestSynth(QuietTestCase):
    """Sanity of the generator itself -- the truth must be trustworthy."""

    def test_deterministic(self):
        np.testing.assert_array_equal(scenario("typical").image, scenario("typical").image)

    def test_star_flux_and_centre(self):
        img = stamp(fwhm=4.0, x=15.3, y=14.6)
        self.assertAlmostEqual(img.sum(), 1e5, delta=1)
        yy, xx = np.indices(img.shape)
        self.assertAlmostEqual((xx * img).sum() / img.sum(), 15.3, places=6)
        self.assertAlmostEqual((yy * img).sum() / img.sum(), 14.6, places=6)

    def test_theta_convention(self):
        # Major axis at 30 deg CCW from +x: second moments must agree.
        img = stamp(fwhm=6.0, ellipticity=0.4, theta=np.deg2rad(30))
        yy, xx = np.indices(img.shape)
        dx, dy = xx - 15, yy - 15
        mxx, myy, mxy = [(img * m).sum() / img.sum() for m in (dx * dx, dy * dy, dx * dy)]
        self.assertAlmostEqual(0.5 * np.arctan2(2 * mxy, mxx - myy), np.deg2rad(30), places=6)

    def test_truth_masks(self):
        f = scenario("typical")
        self.assertEqual(f.truth["mask_hot"].sum(), 3)
        self.assertEqual(f.truth["mask_bad_column"].sum(), 512)
        self.assertTrue(f.truth["mask_line"].any())
        self.assertTrue(f.truth["mask_saturated"].any())
        self.assertLessEqual(f.image.max(), 45000)


class TestStatic(QuietTestCase):
    """Per-stamp measurement helpers on noiseless stars."""

    def test_centroid(self):
        cx, cy = FFS.centroid(stamp(x=15.3, y=14.6))
        self.assertAlmostEqual(cx, 15.3, places=5)
        self.assertAlmostEqual(cy, 14.6, places=5)

    def test_fwhm_round(self):
        fx, fy = FFS.fwhm(stamp(fwhm=4.0))
        self.assertAlmostEqual(fx, 4.0, delta=0.08)
        self.assertAlmostEqual(fy, 4.0, delta=0.08)

    def test_pca_round(self):
        f, e, _ = FFS.pca(stamp(fwhm=4.0))
        self.assertAlmostEqual(f, 4.0, delta=0.05)
        self.assertLess(e, 0.01)

    def test_pca_ellipticity_centred(self):
        _, e, _ = FFS.pca(stamp(fwhm=6.0, ellipticity=0.4, theta=np.deg2rad(30)))
        self.assertAlmostEqual(e, 0.4, delta=0.01)

    def test_pca_ellipticity_off_centre(self):
        # Moments must be taken about the centroid, not the box centre.
        _, e, _ = FFS.pca(stamp(fwhm=6.0, ellipticity=0.3, theta=np.deg2rad(30), x=15.3, y=14.6))
        self.assertAlmostEqual(e, 0.3, delta=0.01)

    def test_pca_theta_convention(self):
        # theta is measured from +x (columns) towards +y (rows).
        for deg in (0, 30, 60, 120, 150):
            with self.subTest(theta=deg):
                _, _, t = FFS.pca(stamp(fwhm=6.0, ellipticity=0.4, theta=np.deg2rad(deg)))
                self.assertLess(axial_diff(t, np.deg2rad(deg)), np.deg2rad(1))

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
        # A skipped (saturated) star must not leave a NaN row or shift the
        # remaining rows.
        stars = [dict(x=100.2, y=100.7, flux=5e6)]
        stars += [dict(x=60 + 60 * i + 0.3, y=300.4, flux=5e4) for i in range(6)]
        f = make_frame(shape=(400, 500), stars=stars, sky=500, saturation=45000, seed=3)
        ffs = FFS(f.image)
        ffs.saturation = 45000
        ffs.calc_frame_fwhm(threshold=10, fwhm=4, box=10, N_stars=50, clip=4)
        self.assertTrue(np.all(np.asarray(ffs.stars["max_adu"]) < 45000))
        self.assertEqual(len(ffs.stars), 6)
        self.assertTrue(np.all(np.isfinite(ffs.stars["fwhm"])))

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
        self.assertAlmostEqual(ffs.max_amplitude, model.max() - model.min(),
                               delta=0.02 * (model.max() - model.min()))

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


if __name__ == "__main__":
    unittest.main()
