import datetime
import unittest

from pyaraucaria.coordinates import *

# >> > parse_deg_or_dms('12:34:56.1')
# 12.58225
# >> > parse_deg_or_dms('12°34′56.1″')
# 12.58225
# >> > parse_deg_or_dms('-12 34 56.1')
# -12.58225
# >> > parse_deg_or_dms('12.58225')


class TestCoo(unittest.TestCase):

    def test_dms_to_float_deg(self):
        tst = [
            ('12°34′56.1″', 12.58225),
            ('12:34:56.1', 12.58225),
            ('-12 34 56.1', -12.58225),
            ('12.58225', 12.58225),
            ('200:30:00', 200.5),
            ]
        res = [deg_to_decimal_deg(d[0]) for d in tst]
        for (_, t), r in zip(tst, res):
            self.assertAlmostEqual(t, r)

    def test_sex_format(self):
        tst = [
            (12.58225, '+12:34:56.1'),
            (-12.58225, '-12:34:56.1'),
            ]
        res = [to_degminsec_sexagesimal(d[0], precision=1) for d in tst]
        for (_, t), r in zip(tst, res):
            self.assertEqual(t, r)


class TestRaDecEpoch(unittest.TestCase):

    def test_ra_dec_epoch(self):
        ra = ra_to_decimal('14:50:42.40')
        dec = dec_to_decimal('74:09:19.73')
        out_ra = ra_to_decimal('14:50:45.27')
        out_dec = dec_to_decimal('74:03:41.46')
        ra_f, dec_f = ra_dec_epoch(ra=ra,
                                   dec=dec,
                                   epoch='2000')
        self.assertAlmostEqual(ra_f, out_ra , places=1)
        self.assertAlmostEqual(dec_f, out_dec, places=1)


class TestRaDecToAzAlt(unittest.TestCase):

    def test_ra_dec_az_alt(self):
        # Altair
        OCA = {'latitude': -24.598056, 'longitude': -70.196389, 'elevation': 2817}
        ra = '19:50:46.68'
        dec = '8:52:03'
        out_az = '304:18:05'
        out_alt = '37:19:12'
        time = datetime.datetime(2023, 3, 27, 15, 0, 16)
        e_az, e_alt = ra_dec_2_az_alt(ra=ra_to_decimal(ra),
                                      dec=dec_to_decimal(dec),
                                      latitude=OCA['latitude'],
                                      longitude=OCA['longitude'],
                                      elevation=OCA['elevation'],
                                      epoch='2000',
                                      time=time)
        out_az = deg_to_decimal_deg(out_az)
        out_alt = deg_to_decimal_deg(out_alt)
        self.assertAlmostEqual(e_az, out_az, places=1)
        self.assertAlmostEqual(e_alt, out_alt, places=1)


class TestSiteSiderealTime(unittest.TestCase):

    def test_site_sidereal_time(self):

        OCA = {'latitude': -24.598056, 'longitude': -70.196389, 'elevation': 2817}
        time = datetime.datetime(2023, 5, 12, 3, 0, 0)
        sidereal_test = '13:37:45.111'
        sidereal_time = site_sidereal_time(longitude=OCA['longitude'],
                                           latitude=OCA['latitude'],
                                           elevation=OCA['elevation'],
                                           time=time)
        sidereal_time_deg = sidereal_time
        sidereal_test_deg = hourangle_to_decimal_deg(sidereal_test)
        self.assertAlmostEqual(sidereal_time_deg, sidereal_test_deg, places=2)


class TestAzAlt2RaDec(unittest.TestCase):

    def test_az_alt_2_ra_dec(self):
        OCA = {'latitude': -24.598056, 'longitude': -70.196389, 'elevation': 2817}
        ra = '19:50:46.68'
        dec = '8:52:03'
        az, alt = ra_dec_2_az_alt(ra=ra_to_decimal(ra),
                                 dec=dec_to_decimal(dec),
                                 latitude=OCA['latitude'],
                                 longitude=OCA['longitude'],
                                 elevation=OCA['elevation'],
                                 epoch='2000')
        ra_2, dec_2 = az_alt_2_ra_dec(az=az,
                                      alt=alt,
                                      latitude=OCA['latitude'],
                                      longitude=OCA['longitude'],
                                      elevation=OCA['elevation'])
        ra_0 = hourangle_to_decimal_deg(ra)
        dec_0 = dec_to_decimal(dec)
        self.assertAlmostEqual(ra_0, ra_2, places=3)
        self.assertAlmostEqual(dec_0, dec_2, places=3)


class TestAzAlt2RaDecAstropy(unittest.TestCase):

    def test_ra_dec_2_az_alt_astropy(self):
        OCA = {'latitude': '-24:35:53', 'longitude': '-70:11:47', 'elevation': '2817'}
        latitude = OCA['latitude']
        longitude = OCA['longitude']
        elevation = OCA['elevation']
        ra = '10:00:00'
        dec = '70:00:00'
        tim = datetime.datetime.now()
        az, alt = ra_dec_2_az_alt(ra=ra_to_decimal(ra), dec=dec_to_decimal(dec), latitude=deg_to_decimal_deg(latitude),
                                  longitude=deg_to_decimal_deg(longitude), elevation=float(elevation), epoch='2000', time=tim)
        ra_2, dec_2 = az_alt_2_ra_dec_astropy(az=az,
                                              alt=alt,
                                              latitude=deg_to_decimal_deg(latitude),
                                              longitude=deg_to_decimal_deg(longitude),
                                              elevation=float(elevation),
                                              epoch='J2000',
                                              calc_time=tim)
        ra_0 = ra_to_decimal(ra)
        dec_0 = dec_to_decimal(dec)
        self.assertAlmostEqual(ra_0, float(ra_2), places=2)
        self.assertAlmostEqual(dec_0, float(dec_2), places=2)


if __name__ == '__main__':
    unittest.main()
