import unittest

import numpy as np

from pyaraucaria.date import *


class TestConvertJd(unittest.TestCase):
    """Tests for converting Julian Date to Barycentric Julian Date."""

    def test_to_bjd_float(self):
        """Convert a fixed OCA observation epoch and sky position to expected BJD."""
        oca_loc = {'latitude': -24.59806, 'longitude': -70.19638, 'elevation': 2817.0}
        jd = 2460000
        ra = 10.0
        dec = 70.0
        expected_bjd = 2460000.000057626819

        bjd = jd_to_bjd(
            jd=jd,
            obj_ra=ra,
            obj_dec=dec,
            observ_lat=oca_loc['latitude'],
            observ_lon=oca_loc['longitude'],
            observ_elev=oca_loc['elevation']
        )


        self.assertIsInstance(bjd, float)
        self.assertAlmostEqual(bjd, expected_bjd, places=8)

    def test_to_bjd_array(self):
        """Convert a fixed OCA observation epoch and sky position to expected BJD."""
        oca_loc = {'latitude': -24.59806, 'longitude': -70.19638, 'elevation': 2817.0}
        jd = np.array([
            2450000,
            2455000,
            2460000
        ])
        ra = 10.0
        dec = 70.0
        expected_bjd = np.array([
            2450000.003263008708,
            2454999.998134490116,
            2460000.000057626819

            ])

        bjd = jd_to_bjd(
            jd=jd,
            obj_ra=ra,
            obj_dec=dec,
            observ_lat=oca_loc['latitude'],
            observ_lon=oca_loc['longitude'],
            observ_elev=oca_loc['elevation']
        )

        # self.assertIsInstance(bjd, np.ndarray)
        # self.assertAlmostEqual(bjd, expected_bjd, places=8)

        for n, m in enumerate(expected_bjd):
            self.assertAlmostEqual(bjd[n], expected_bjd[n], places=8)


class TestHmsToDays(unittest.TestCase):

    def test_hms_to_days(self):
        self.assertEqual(hms_to_days(), 0)
        self.assertAlmostEqual(hms_to_days(hours=12), 0.5)
        self.assertAlmostEqual(hms_to_days(minutes=30), 1 / 48)
        self.assertAlmostEqual(hms_to_days(seconds=60), 1 / 1440)
        self.assertAlmostEqual(hms_to_days(hours=1, minutes=2, seconds=3), 3723 / 86400)


if __name__ == '__main__':
    unittest.main()
