"""Tests for FFS.mk_pixels_map / mk_hot_pixels_mask / mk_cold_pixels_mask (calibration frames)."""

import unittest

import numpy as np

from pyaraucaria.ffs import FFS

N = 256


class TestHotColdPixels(unittest.TestCase):

    def setUp(self):
        self.rng = np.random.default_rng(1)
        self.yy, self.xx = np.mgrid[0:N, 0:N]

    def test_dark_hot_pixels(self):
        img = 1000 + self.rng.normal(0, 5, (N, N))
        hot = self.rng.integers(0, N, (50, 2))
        img[hot[:, 0], hot[:, 1]] += 200
        ffs = FFS(img)
        ffs.mk_pixels_map()
        ffs.mk_hot_pixels_mask(threshold=50)
        ffs.mk_cold_pixels_mask(threshold=50)
        truth = np.zeros((N, N), bool)
        truth[hot[:, 0], hot[:, 1]] = True
        self.assertTrue(np.array_equal(ffs.masks["hot_pixels"], truth))
        self.assertFalse(ffs.masks["cold_pixels"].any())
        self.assertIn("hot_pixels", ffs.exclude)
        self.assertIn("cold_pixels", ffs.exclude)

    def test_flat_cold_pixels_with_vignetting(self):
        vign = 30000 * (1 - 0.3 * (np.hypot(self.xx - N / 2, self.yy - N / 2) / (N / 2)) ** 2)
        sens = np.ones((N, N))
        cold = self.rng.integers(0, N, (50, 2))
        sens[cold[:, 0], cold[:, 1]] = 0.8
        ffs = FFS(self.rng.poisson(vign * sens).astype(float))
        ffs.mk_pixels_map()
        ffs.mk_cold_pixels_mask(threshold=2000)
        ffs.mk_hot_pixels_mask(threshold=2000)
        truth = np.zeros((N, N), bool)
        truth[cold[:, 0], cold[:, 1]] = True
        self.assertTrue(np.array_equal(ffs.masks["cold_pixels"], truth))
        self.assertFalse(ffs.masks["hot_pixels"].any())

    def test_map_in_adu(self):
        img = np.full((N, N), 1000.0)
        img[50, 60] = 1300
        img[70, 80] = 900
        ffs = FFS(img)
        ffs.mk_pixels_map()
        self.assertEqual(ffs.maps["pixels"][50, 60], 300)
        self.assertEqual(ffs.maps["pixels"][70, 80], -100)


if __name__ == "__main__":
    unittest.main()
