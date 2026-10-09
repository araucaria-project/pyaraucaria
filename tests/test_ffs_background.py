"""find_stars / star_info with the sky map and exclusion masks (bad_mask)."""

import unittest

import numpy as np

from pyaraucaria.ffs import FFS
from tests.ffs_synth import make_frame, star_field


def matched(coo, truth, tol=2.0):
    if len(coo) == 0:
        return 0
    d = np.hypot(*(coo[:, None, :] - truth[None, :, :]).transpose(2, 0, 1)).min(axis=1)
    return int((d < tol).sum())


class TestFindStarsBackground(unittest.TestCase):

    def setUp(self):
        self.stars = star_field(n=60, seed=5)
        self.truth = np.array([[s["y"], s["x"]] for s in self.stars])
        f = make_frame(stars=self.stars, sky=1000, seed=6)
        xx = np.mgrid[0:f.image.shape[0], 0:f.image.shape[1]][1]
        self.image = f.image + 3.0 * xx        # ~1500 ADU gradient across the frame

    def test_sky_map_recovers_stars_on_gradient(self):
        ffs = FFS(self.image)
        ffs.mk_stats()
        ffs.find_stars(threshold=5, fwhm=4)
        n_median = matched(ffs.coo, self.truth)

        ffs = FFS(self.image)
        ffs.mk_stats()
        ffs.sky_map()
        ffs.find_stars(threshold=5, fwhm=4)
        n_sky = matched(ffs.coo, self.truth)

        self.assertGreater(n_sky, n_median)
        self.assertGreaterEqual(n_sky, 45)
        self.assertLessEqual(len(ffs.coo) - n_sky, 2)

    def test_bad_pixels_are_not_stars(self):
        img = make_frame(stars=self.stars, sky=1000, seed=6).image
        hot = [(100, 37), (300, 211), (450, 480)]
        for y, x in hot:
            img[y, x] += 20000
        ffs = FFS(img)
        ffs.masks["hot"] = np.zeros(img.shape, bool)
        for y, x in hot:
            ffs.masks["hot"][y, x] = True
        ffs.exclude.add("hot")
        ffs.mk_stats()
        ffs.find_stars(threshold=5, fwhm=4)
        found = {tuple(c) for c in ffs.coo.tolist()}
        self.assertFalse(found & set(hot))

    def test_nan_pixels_do_not_break_find_stars(self):
        img = make_frame(stars=self.stars, sky=1000, seed=6).image
        img[:, :20] = np.nan
        ffs = FFS(img)
        ffs.mk_stats()
        ffs.find_stars(threshold=5, fwhm=4)
        self.assertGreaterEqual(matched(ffs.coo, self.truth), 45)
        self.assertTrue(np.all(np.isfinite(ffs.adu)))

    def test_star_info_skips_stars_with_bad_pixels(self):
        img = make_frame(stars=self.stars, sky=1000, seed=6).image
        ffs = FFS(img)
        ffs.mk_stats()
        ffs.find_stars(threshold=5, fwhm=4)
        y, x = ffs.coo[0]
        ffs.masks["bad"] = np.zeros(img.shape, bool)
        ffs.masks["bad"][y, x + 3] = True           # one bad pixel next to the brightest star
        ffs.exclude.add("bad")
        n_before = len(ffs.coo)
        ffs.star_info(box=8)
        # rejected stars are dropped from coo; the others are measured
        self.assertNotIn([y, x], ffs.coo.tolist())
        self.assertGreater(len(ffs.coo), n_before // 2)
        self.assertTrue(np.all(np.isfinite(ffs.fwhm)))


if __name__ == "__main__":
    unittest.main()
