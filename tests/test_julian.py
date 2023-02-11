import unittest
from pyaraucaria.date import *


class TestJulian(unittest.TestCase):

    def test_to_julian(self):
        dt = [
            ('2010.10.02 4:27:15.261', 2455471.685593298),
            ('2010-10-02T04:27:15.261', 2455471.685593298),
            ('31/01/2022', 2459610.5),
            ('01/02/2022', 2459611.5),
            ]
        jd1 = [datetime_to_julian(d[0]) for d in dt]
        # jd2 = [datetime_to_julian_old(d[0]) for d in dt]
        for i, j in enumerate(jd1):
            self.assertAlmostEqual(j, dt[i][1])

    def test_julian_to_iso(self):
        self.assertEqual(julian_to_iso(2455471.685593298), "2010-10-02T04:27:15.261")
        self.assertEqual(julian_to_iso(2459499.166666667), "2021-10-11T16:00:00.000")
        self.assertEqual(julian_to_iso(2459499.0),         "2021-10-11T12:00:00.000")
        self.assertEqual(julian_to_iso(2459499.5),         "2021-10-12T00:00:00.000")

    def test_fromat_iso(self):
        self.assertEqual(tuple_to_iso(2021, 5, 12, 1, 12, 3.9191, seconds_decimals=0), "2021-05-12T01:12:04")
        self.assertEqual(tuple_to_iso(2021, 5, 12, 1, 12, 3.9191, seconds_decimals=2), "2021-05-12T01:12:03.92")
        self.assertEqual(tuple_to_iso(2021, 5, 12, 1, 12, 3.9191, seconds_decimals=4), "2021-05-12T01:12:03.9191")
        self.assertEqual(tuple_to_iso(2021, 5, 12, 1, 12, 3.9191, include_time=False), "2021-05-12")
        self.assertEqual(tuple_to_iso(2021, 5, 12, 1, 12, 3.9191, include_date=False), "01:12:03.919")
        self.assertEqual(tuple_to_iso(2021, 5, 12, 1, 12, 3.9191, separator=' '), "2021-05-12 01:12:03.919")

    def test_reduced_julian(self):
        self.assertEqual(julian_to_iso(2455471.685593298), "2010-10-02T04:27:15.261")
        self.assertEqual(julian_to_iso(  55471.685593298), "2010-10-02T04:27:15.261")
        self.assertEqual(julian_to_iso(   5471.685593298), "2010-10-02T04:27:15.261")

    @unittest.skip("Not implemented")
    def test_hjd(self):
        jd = datetime_to_julian("2010-10-02T04:27:15.261")
        ra, dec = "", ""
        hcor = helio_corr(jd, ra, dec)
        hjd = jd + hcor

if __name__ == '__main__':
    unittest.main()
