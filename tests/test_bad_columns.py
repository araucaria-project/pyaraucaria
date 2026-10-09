"""Tests for FFS.mk_columns_map / mk_bad_columns_mask on synthetic calibration frames (no stars)."""

import unittest

import numpy as np

from pyaraucaria.ffs import FFS

N = 512


def found(ffs):
    cols = sorted(set(np.nonzero(ffs.masks["bad_columns"])[1].tolist()))
    rows = sorted(set(np.nonzero(ffs.masks["bad_rows"])[0].tolist()))
    return cols, rows


class TestBadColumns(unittest.TestCase):

    def setUp(self):
        self.rng = np.random.default_rng(1)
        self.yy, self.xx = np.mgrid[0:N, 0:N]

    def test_maps_in_adu(self):
        img = 1000 + self.rng.normal(0, 5, (N, N))
        img[:, 100] += 10
        ffs = FFS(img)
        ffs.mk_columns_map()
        self.assertEqual(ffs.maps["columns"].shape, img.shape)
        self.assertAlmostEqual(np.median(ffs.maps["columns"][:, 100]), 10, delta=1.5)

    def test_clean_bias_has_no_bad_columns(self):
        ffs = FFS(np.round(1000 + self.rng.normal(0, 5, (N, N))).astype(np.uint16))
        ffs.mk_columns_map()
        ffs.mk_bad_columns_mask(threshold=5)
        self.assertEqual(found(ffs), ([], []))
        self.assertIn("bad_columns", ffs.exclude)
        self.assertIn("bad_rows", ffs.exclude)

    def test_bias_columns_and_row(self):
        img = 1000 + self.rng.normal(0, 5, (N, N))
        img[:, 100] += 10
        img[:, 300:302] -= 10       # 2 px wide, dark
        img[200, :] += 10
        ffs = FFS(np.round(img).astype(np.uint16))
        ffs.mk_columns_map()
        ffs.mk_bad_columns_mask(threshold=5)
        self.assertEqual(found(ffs), ([100, 300, 301], [200]))

    def test_partial_column_in_dark_with_cosmics(self):
        img = 1000 + self.rng.normal(0, 5, (N, N))
        hits = self.rng.integers(0, N, (500, 2))
        img[hits[:, 0], hits[:, 1]] += self.rng.uniform(100, 5000, 500)
        img[256:, 150] += 20        # defect from row 256 down
        ffs = FFS(img)
        ffs.mk_columns_map()
        ffs.mk_bad_columns_mask(threshold=5)
        self.assertEqual(found(ffs), ([150], []))
        self.assertTrue(ffs.masks["bad_columns"][256:, 150].all())
        self.assertFalse(ffs.masks["bad_columns"][:256, 150].any())

    def test_flat_with_vignetting(self):
        vign = 30000 * (1 - 0.3 * (np.hypot(self.xx - N / 2, self.yy - N / 2) / (N / 2)) ** 2)
        sens = np.ones((N, N))
        sens[:, 100] *= 0.98
        sens[:, 400:402] *= 0.98
        sens[300, :] *= 0.98
        ffs = FFS(self.rng.poisson(vign * sens).astype(float))
        ffs.sky_map(n_segments=32)  # vignetting needs a finer sky map than the default
        ffs.mk_columns_map()
        ffs.mk_bad_columns_mask(threshold=300)
        self.assertEqual(found(ffs), ([100, 400, 401], [300]))


if __name__ == "__main__":
    unittest.main()
