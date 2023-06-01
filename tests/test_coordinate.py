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

    def test_ra_dec_epoch_1(self):
        RUD = {'latitude': '51:17:05', 'longitude': '22:35:50', 'elevation': '235'}
        ra = '5:14:32.30'
        dec = '-8:12:06'
        out_ra = '5:15:38.46'
        out_dec = '-8:10:38.83'
        epoch = ra_dec_epoch(ra=ra, dec=dec,
                             latitude=RUD['latitude'], longitude=RUD['longitude'],
                             elevation=RUD['elevation'], epoch='2000')
        e_ra = ra_to_decimal(out_ra)
        e_dec = dec_to_decimal(out_dec)
        out_epoch_ra = ra_to_decimal(epoch[0])
        out_epoch_dec = dec_to_decimal(epoch[1])
        self.assertAlmostEqual(e_ra, out_epoch_ra, places=2)
        self.assertAlmostEqual(e_dec, out_epoch_dec, places=2)

class TestRaDecToAzAlt(unittest.TestCase):

    def test_ra_dec_az_alt(self):
        # Altair
        OCA = {'latitude': '-24:35:53', 'longitude': '-70:11:47', 'elevation': '2817'}
        ra = '19:50:46.68'
        dec = '8:52:03'
        out_az = '304:18:05'
        out_alt = '37:19:12'
        time = datetime.datetime(2023, 3, 27, 15, 0, 16)
        az_alt = ra_dec_2_az_alt(ra=ra, dec=dec,
                                 latitude=OCA['latitude'], longitude=OCA['longitude'],
                                 elevation=OCA['elevation'], epoch='2000', time=time)


        e_az = deg_str_to_deg(out_az)
        e_alt = deg_str_to_deg(out_alt)
        out_az = deg_str_to_deg(az_alt[0])
        out_alt = deg_str_to_deg(az_alt[1])
        self.assertAlmostEqual(e_az, out_az, places=1)
        self.assertAlmostEqual(e_alt, out_alt, places=1)


class TestStrToDeg(unittest.TestCase):

    def test_str_to_deg(self):
        a = '10:30:00'
        b = 10.5
        a1 = '-05:00:30.30'
        b1 = -5.008416667
        aa = deg_str_to_deg(a)
        aa1 = deg_str_to_deg(a1)
        self.assertEqual(aa, b)
        self.assertAlmostEqual(aa1, b1, 8)


class TestSiteSiderealTime(unittest.TestCase):

    def test_site_sidereal_time(self):

        OCA = {'latitude': '-24:35:53', 'longitude': '-70:11:47', 'elevation': '2817'}
        time = datetime.datetime(2023, 5, 12, 3, 0, 0)
        sidereal_test = '13:37:45.111'
        sidereal_time = site_sidereal_time(longitude=OCA['longitude'], latitude=OCA['latitude'],
                    elevation=OCA['elevation'], time=time)
        sidereal_time_deg_2 = site_sidereal_time(longitude=OCA['longitude'], latitude=OCA['latitude'],
                                           elevation=OCA['elevation'], time=time, deg_output=True)
        sidereal_time_deg = hourangle_to_decimal_deg(sidereal_time)
        sidereal_test_deg = hourangle_to_decimal_deg(sidereal_test)
        self.assertAlmostEqual(sidereal_time_deg, sidereal_test_deg, places=2)
        self.assertAlmostEqual(sidereal_time_deg_2, sidereal_test_deg, places=2)

class TestAzAlt2RaDec(unittest.TestCase):

    def test_az_alt_2_ra_dec(self):
        OCA = {'latitude': '-24:35:53', 'longitude': '-70:11:47', 'elevation': '2817'}
        latitude = OCA['latitude']
        longitude = OCA['longitude']
        elevation = OCA['elevation']
        ra = '19:50:46.68'
        dec = '8:52:03'
        az, alt = ra_dec_2_az_alt(ra=ra, dec=dec, latitude=latitude,
                                  longitude=longitude, elevation=elevation, epoch='2000')
        ra_2, dec_2 = az_alt_2_ra_dec(deg_str_to_deg(az), deg_str_to_deg(alt), latitude=deg_str_to_deg(OCA['latitude']),
                                      longitude=deg_str_to_deg(OCA['longitude']), elevation=float(OCA['elevation']))
        ra_0 = hourangle_to_decimal_deg(ra)
        dec_0 = dec_to_decimal(dec)
        self.assertAlmostEqual(ra_0, ra_2, places=3)
        self.assertAlmostEqual(dec_0, dec_2, places=3)

if __name__ == '__main__':
    unittest.main()
