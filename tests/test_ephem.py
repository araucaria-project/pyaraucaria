import unittest
from datetime import datetime

from astropy.coordinates import EarthLocation

from pyaraucaria.ephem import oca_location, as_location, get_moon_ephem


class TestLocation(unittest.TestCase):

    def test_as_location(self):
        x = [
            oca_location,
            EarthLocation.of_site('OCA'),
            'oca',
            EarthLocation.from_geodetic(lon=-70.20128, lat=-24.59867, height=2817),
            (-70.20128, -24.59867, 2817),
        ]
        for loc in x:
            self.assertEqual(as_location(loc), oca_location)

class TestEphem(unittest.TestCase):

    def test_moon(self):
        e2 = get_moon_ephem([datetime(2021, 1, 1), datetime(2021, 1, 2)])
        self.assertEqual(len(e2), 2)
        e = get_moon_ephem()
        self.assertEqual(len(e), 1)
        self.assertGreater(e[0]['alt'], -90)


if __name__ == '__main__':
    unittest.main()
