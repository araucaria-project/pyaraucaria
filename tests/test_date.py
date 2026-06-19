import unittest

import astropy.units as u
from astropy.coordinates import EarthLocation, SkyCoord

from pyaraucaria.date import jd_to_bjd


class TestJdToBjd(unittest.TestCase):

    def test_jd_to_bjd(self):
        oca_loc = {'latitude': -24.59806, 'longitude': -70.19638, 'elevation': 2817.0}
        bjd_result = 2460000.000057627
        jd = 2460000
        ra = 10.0
        dec = 70.0
        bjd = jd_to_bjd(
            jd=jd,
            obj_ra=ra,
            obj_dec=dec,
            observ_lat=oca_loc['latitude'],
            observ_lon=oca_loc['longitude'],
            observ_elev=oca_loc['elevation']

        )
        self.assertAlmostEqual(bjd, bjd_result, places=8)

    def test_jd_to_bjd_with_skycoord_and_earthlocation(self):
        jd = 2460000
        target = SkyCoord(ra=10.0 * u.deg, dec=70.0 * u.deg)
        location = EarthLocation(lat=-24.59806 * u.deg, lon=-70.19638 * u.deg, height=2817.0 * u.m)

        bjd_numeric = jd_to_bjd(jd, 10.0, 70.0, -24.59806, -70.19638, 2817.0)
        bjd_objects = jd_to_bjd(jd=jd, obj_ra=target, obj_dec=None, observ_lat=location)

        self.assertAlmostEqual(bjd_objects, bjd_numeric, places=12)


if __name__ == '__main__':
    unittest.main()
