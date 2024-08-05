import datetime
import unittest
from pyaraucaria.ephemeris import moon_phase


class TestAzAlt2RaDecAstropy(unittest.TestCase):

    def test_ra_dec_2_az_alt_astropy(self):
        OCA = {'latitude': -24.59855, 'longitude': -70.20126, 'elevation': 2817}
        latitude = OCA['latitude']
        longitude = OCA['longitude']
        elevation = OCA['elevation']
        tim = datetime.datetime.fromisoformat('2024-08-11T00:00:00')
        test_phase = 34.6
        mp = moon_phase(date_utc=tim, latitude=latitude, longitude=longitude, elevation=elevation)
        self.assertAlmostEqual(mp, test_phase, places=1)


if __name__ == '__main__':
    unittest.main()