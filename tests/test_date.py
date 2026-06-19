import datetime
import unittest

from pyaraucaria.date import *


class TestJdToBjd(unittest.TestCase):

    def test_jd_to_bjd(self):
        oca_loc = {'latitude': -24.59806, 'longitude': -70.19638, 'elevation': 2817.0}
        latitude = oca_loc['latitude']
        longitude = oca_loc['longitude']
        elevation = oca_loc['elevation']
        bjd_result = 2460000.000057627
        jd = 2460000
        ra = 10.0
        dec = 70.0
        tim = datetime.datetime.now()
        bjd = jd_to_bjd(
            jd=jd,
            obj_ra=ra,
            obj_dec=dec,
            observ_lat=oca_loc['latitude'],
            observ_lon=oca_loc['longitude'],
            observ_elev=oca_loc['elevation']

        )
        self.assertAlmostEqual(bjd, bjd_result, places=8)


if __name__ == '__main__':
    unittest.main()
