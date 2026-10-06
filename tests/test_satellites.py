"""Truth tests for satellite-trail detection (``pyaraucaria.satellites``).

Every scenario injects known trails (``tests/ffs_synth.py``) and checks that
each trail is found once, with the right geometry, that its pixel mask covers
the truth mask, and that frames without trails give no detections.
"""

import os
import time
import unittest
import warnings

import numpy as np
from scipy.ndimage import binary_dilation

from pyaraucaria.ffs import FFS
from pyaraucaria.satellites import find_satellites
from tests.ffs_synth import make_frame, scenario, star_field

JK15C = os.path.expanduser("~/work/data/fits_examples/jk15c_0154_73367.fits")


def trail(x0, y0, x1, y1, flux, width=2.5):
    return dict(x0=x0, y0=y0, x1=x1, y1=y1, flux=flux, width=width)


def line_params(ln):
    """(theta, rho) of the infinite line through a truth segment, FFS convention."""
    dx, dy = ln["x1"] - ln["x0"], ln["y1"] - ln["y0"]
    n = np.array([-dy, dx]) / np.hypot(dx, dy)
    theta = np.arctan2(n[1], n[0])
    rho = ln["x0"] * n[0] + ln["y0"] * n[1]
    # fold into theta in [-pi/2, pi/2)
    if theta >= np.pi / 2:
        theta, rho = theta - np.pi, -rho
    elif theta < -np.pi / 2:
        theta, rho = theta + np.pi, -rho
    return theta, rho


def axial_diff(a, b):
    return np.abs((a - b + np.pi / 2) % np.pi - np.pi / 2)


def visible_segment(ln, shape):
    """Truth segment clipped to the frame."""
    ny, nx = shape
    t = np.linspace(0, 1, 4001)
    x = ln["x0"] + t * (ln["x1"] - ln["x0"])
    y = ln["y0"] + t * (ln["y1"] - ln["y0"])
    ok = (x >= 0) & (x <= nx - 1) & (y >= 0) & (y <= ny - 1)
    return x[ok][0], y[ok][0], x[ok][-1], y[ok][-1]


def match_trail(table, ln, shape, max_angle=np.deg2rad(1.0), max_offset=3.0):
    """Rows of ``table`` (satellites only) that describe truth trail ``ln``."""
    theta, _ = line_params(ln)
    xm, ym = (ln["x0"] + ln["x1"]) / 2, (ln["y0"] + ln["y1"]) / 2
    rows = []
    for i, row in enumerate(table):
        if row["kind"] != "satellite":
            continue
        off = abs(xm * np.cos(row["theta"]) + ym * np.sin(row["theta"]) - row["rho"])
        if axial_diff(row["theta"], theta) < max_angle and off < max_offset:
            rows.append(i)
    return rows


def endpoint_error(row, ln, shape):
    xa, ya, xb, yb = visible_segment(ln, shape)
    d1 = max(np.hypot(row["x0"] - xa, row["y0"] - ya), np.hypot(row["x1"] - xb, row["y1"] - yb))
    d2 = max(np.hypot(row["x0"] - xb, row["y0"] - yb), np.hypot(row["x1"] - xa, row["y1"] - ya))
    return min(d1, d2)


class SatelliteTestCase(unittest.TestCase):

    def setUp(self):
        self._w = warnings.catch_warnings()
        self._w.__enter__()
        warnings.simplefilter("ignore", RuntimeWarning)

    def tearDown(self):
        self._w.__exit__(None, None, None)

    def detect(self, frame, **kw):
        return find_satellites(frame.image, **kw)

    def assert_trails_found(self, frame, table, mask, endpoint_tol=12.0,
                            completeness=0.95, purity=0.90):
        sats = [r for r in table if r["kind"] == "satellite"]
        truth_lines = frame.truth["lines"]
        self.assertEqual(len(sats), len(truth_lines),
                         f"expected {len(truth_lines)} satellites, got {len(sats)}:\n{table}")
        for ln in truth_lines:
            rows = match_trail(table, ln, frame.image.shape)
            self.assertEqual(len(rows), 1, f"trail {ln} matched rows {rows}:\n{table}")
            self.assertLess(endpoint_error(table[rows[0]], ln, frame.image.shape), endpoint_tol,
                            f"endpoints off for {ln}:\n{table[rows[0]]}")
        truth = frame.truth["mask_line"]
        covered = (mask & truth).sum() / truth.sum()
        self.assertGreaterEqual(covered, completeness, f"mask covers {covered:.3f} of truth")
        near = binary_dilation(truth, iterations=3)
        pure = (mask & near).sum() / max(mask.sum(), 1)
        self.assertGreaterEqual(pure, purity, f"only {pure:.3f} of mask is near a trail")


class TestDetection(SatelliteTestCase):

    def test_bright_trail_with_stars(self):
        f = make_frame(stars=star_field(n=40, seed=11), sky=1000,
                       lines=[trail(0, 60, 511, 300, 300)], seed=12)
        self.assert_trails_found(f, *self.detect(f))

    def test_faint_trail(self):
        # 1.25 sigma per pixel at the trail centre
        f = make_frame(stars=star_field(n=40, seed=13), sky=1000,
                       lines=[trail(30, 0, 480, 511, 40, width=3.0)], seed=14)
        self.assert_trails_found(f, *self.detect(f))

    def test_short_segment(self):
        f = make_frame(stars=star_field(n=30, seed=15), sky=1000,
                       lines=[trail(100, 120, 330, 210, 150)], seed=16)
        self.assert_trails_found(f, *self.detect(f), endpoint_tol=10.0)

    def test_two_crossing_trails(self):
        f = make_frame(stars=star_field(n=30, seed=17), sky=1000,
                       lines=[trail(0, 50, 511, 450, 200), trail(0, 400, 511, 120, 120)], seed=18)
        self.assert_trails_found(f, *self.detect(f))

    def test_trail_through_saturated_star(self):
        stars = star_field(n=30, seed=19) + [dict(x=256.4, y=256.2, flux=4e6, fwhm=4.0)]
        f = make_frame(stars=stars, sky=1000, saturation=45000,
                       lines=[trail(0, 130, 511, 383, 150)], seed=20)
        table, mask = self.detect(f)
        self.assert_trails_found(f, table, mask, purity=0.85)

    def test_near_vertical_trail_is_satellite(self):
        # 3 deg from vertical and 3 px wide: a satellite, not a bad column.
        f = make_frame(stars=star_field(n=30, seed=21), sky=1000,
                       lines=[trail(200, 0, 226.8, 511, 200, width=3.0)], seed=22)
        self.assert_trails_found(f, *self.detect(f))

    def test_satellite_train_stays_satellites(self):
        # Several bright parallel trails (e.g. a Starlink train) must not be
        # taken for a faint periodic pattern.
        lines = [trail(0, 40 + 60 * i, 511, 200 + 60 * i, 250, width=2.5) for i in range(4)]
        f = make_frame(stars=star_field(n=30, seed=31), sky=1000, lines=lines, seed=32)
        self.assert_trails_found(f, *self.detect(f))

    def test_typical_scenario(self):
        # One trail plus a hot column, hot pixels and saturated stars.
        f = scenario("typical")
        table, mask = self.detect(f)
        self.assert_trails_found(f, table, mask)
        # column x=260 may be masked only where the trail crosses it
        near_trail = binary_dilation(f.truth["mask_line"], iterations=10)[:, 260]
        self.assertFalse(mask[:, 260][~near_trail].any(),
                         "bad column must not be masked as a satellite")
        self.assertEqual(list(table["kind"]).count("column"), 1)

    def test_larger_frame_binned(self):
        shape = (2048, 2048)
        f = make_frame(shape=shape, stars=star_field(n=300, shape=shape, seed=23), sky=1000,
                       lines=[trail(0, 300, 2047, 1700, 30, width=3.0)], seed=24)
        t = time.time()
        table, mask = self.detect(f)
        elapsed = time.time() - t
        self.assert_trails_found(f, table, mask, endpoint_tol=20.0)
        self.assertLess(elapsed, 10.0)


class TestNoFalsePositives(SatelliteTestCase):

    def test_scenarios_without_trails(self):
        for name in ("clean", "elliptical", "sparse", "empty", "gradient"):
            with self.subTest(scenario=name):
                f = scenario(name)
                table, mask = self.detect(f)
                self.assertEqual(sum(r["kind"] == "satellite" for r in table), 0, f"\n{table}")
                self.assertFalse(mask.any())

    def test_bad_column_is_not_a_satellite(self):
        f = make_frame(stars=star_field(n=30, seed=25), sky=1000,
                       bad_columns=[dict(x=200, kind="hot", value=300.0)], seed=26)
        table, mask = self.detect(f)
        self.assertEqual(sum(r["kind"] == "satellite" for r in table), 0, f"\n{table}")
        self.assertFalse(mask.any())

    def test_faint_periodic_pattern(self):
        # Faint parallel bands across the frame (0.94 sigma/px, every 100 px),
        # like the flat/detector pattern seen on jk15c: detected, but
        # reported as "pattern" and never masked.
        shape = (1024, 1024)
        bands = [trail(-200 + 100 * i, 0, -200 + 100 * i + 170, 1023, 30.0, width=4.0) for i in range(12)]
        f = make_frame(shape=shape, stars=star_field(n=80, shape=shape, seed=33), sky=1000,
                       lines=bands, seed=34)
        table, mask = self.detect(f)
        self.assertGreaterEqual(sum(r["kind"] == "pattern" for r in table), 3, f"\n{table}")
        self.assertEqual(sum(r["kind"] == "satellite" for r in table), 0, f"\n{table}")
        self.assertFalse(mask.any())

    def test_saturated_stars_with_bleed(self):
        stars = star_field(n=30, seed=27) + [dict(x=150.3, y=380.7, flux=5e6, fwhm=4.0),
                                             dict(x=400.6, y=120.2, flux=4e6, fwhm=4.0)]
        f = make_frame(stars=stars, sky=1000, saturation=45000, seed=28)
        table, mask = self.detect(f)
        self.assertEqual(sum(r["kind"] == "satellite" for r in table), 0, f"\n{table}")

    @unittest.skipUnless(os.path.exists(JK15C), "OCA example frame not available")
    def test_real_frame_without_trails(self):
        from astropy.io import fits
        data = fits.getdata(JK15C)
        t = time.time()
        table, mask = find_satellites(data)
        elapsed = time.time() - t
        self.assertEqual(sum(r["kind"] == "satellite" for r in table), 0, f"\n{table}")
        self.assertLess(elapsed, 15.0)


class TestFFSIntegration(SatelliteTestCase):

    def test_ffs_method_stores_results(self):
        f = scenario("typical")
        ffs = FFS(f.image)
        ffs.find_satellites()
        self.assertEqual(ffs.satellite_mask.shape, f.image.shape)
        self.assertEqual(int(np.sum(ffs.satellites["kind"] == "satellite")), 1)
        self.assertIn("satellites", ffs.stats)
        self.assertEqual(set(ffs.stats["satellites"]), set(ffs.stats_description["satellites"]))


if __name__ == "__main__":
    unittest.main()
