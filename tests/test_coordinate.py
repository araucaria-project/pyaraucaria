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

    def test_ra_dec_epoch_precesses(self):
        """Precessing from an older epoch (J2010) changes coords measurably vs J2000 input."""
        ra = ra_to_decimal('19:50:46.68')  # Altair
        dec = dec_to_decimal('8:52:03')
        ra_from_j2000, dec_from_j2000 = ra_dec_epoch(ra=ra, dec=dec, epoch='2000')
        ra_from_j2010, dec_from_j2010 = ra_dec_epoch(ra=ra, dec=dec, epoch='2010')
        # The same input RA/Dec at two different epochs should yield different "now" positions
        self.assertNotAlmostEqual(ra_from_j2000, ra_from_j2010, places=3)


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


class TestRaDecToAzAltPrecision(unittest.TestCase):
    """Test ra_dec_2_az_alt against astropy-computed reference values at arcsecond level."""

    OCA = {'latitude': -24.598056, 'longitude': -70.196389, 'elevation': 2817}

    def test_altair_precision(self):
        """Altair at OCA 2023-03-27 15:00:16 — reference from FK5(J2000)→AltAz transform.

        astropy reference: az=304.2969°, alt=37.3196°
        """
        time = datetime.datetime(2023, 3, 27, 15, 0, 16)
        az, alt = ra_dec_2_az_alt(
            ra=ra_to_decimal('19:50:46.68'),
            dec=dec_to_decimal('8:52:03'),
            latitude=self.OCA['latitude'],
            longitude=self.OCA['longitude'],
            elevation=self.OCA['elevation'],
            epoch='2000',
            time=time,
        )
        self.assertAlmostEqual(az, 304.2969, places=2)
        self.assertAlmostEqual(alt, 37.3196, places=2)

    def test_canopus_az_alt(self):
        """Canopus (alpha Car) at OCA 2023-03-27 15:00:16 — independent accuracy check.

        J2000 coords: RA=06:23:57.11, Dec=-52:41:44.4
        astropy reference: az=147.019°, alt=4.918°
        """
        time = datetime.datetime(2023, 3, 27, 15, 0, 16)
        az, alt = ra_dec_2_az_alt(
            ra=ra_to_decimal('6:23:57.11'),
            dec=dec_to_decimal('-52:41:44.4'),
            latitude=self.OCA['latitude'],
            longitude=self.OCA['longitude'],
            elevation=self.OCA['elevation'],
            epoch='2000',
            time=time,
        )
        self.assertAlmostEqual(az, 147.019, places=2)
        self.assertAlmostEqual(alt, 4.918, places=2)


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

    def test_sidereal_time_precision(self):
        """OCA 2023-05-12 03:00:00 — reference from astropy apparent LAST.

        astropy LAST: 204.43521° = 13h 37m 44.45s
        """
        OCA = {'longitude': -70.196389, 'latitude': -24.598056, 'elevation': 2817}
        time = datetime.datetime(2023, 5, 12, 3, 0, 0)
        lst = site_sidereal_time(
            longitude=OCA['longitude'],
            latitude=OCA['latitude'],
            elevation=OCA['elevation'],
            time=time,
        )
        self.assertAlmostEqual(lst, 204.43521, places=3)


class TestAzAlt2RaDec(unittest.TestCase):

    def test_az_alt_2_ra_dec(self):
        OCA = {'latitude': -24.598056, 'longitude': -70.196389, 'elevation': 2817}
        ra = '19:50:46.68'
        dec = '8:52:03'
        fixed_time = datetime.datetime(2023, 3, 27, 15, 0, 16)
        az, alt = ra_dec_2_az_alt(ra=ra_to_decimal(ra),
                                 dec=dec_to_decimal(dec),
                                 latitude=OCA['latitude'],
                                 longitude=OCA['longitude'],
                                 elevation=OCA['elevation'],
                                 epoch='2000',
                                 time=fixed_time)
        ra_2, dec_2 = az_alt_2_ra_dec(az=az,
                                      alt=alt,
                                      latitude=OCA['latitude'],
                                      longitude=OCA['longitude'],
                                      elevation=OCA['elevation'],
                                      time=fixed_time)
        ra_0 = hourangle_to_decimal_deg(ra)
        dec_0 = dec_to_decimal(dec)
        self.assertAlmostEqual(ra_0, ra_2, places=3)
        self.assertAlmostEqual(dec_0, dec_2, places=3)


class TestAzAlt2RaDecAccuracy(unittest.TestCase):
    """Test az_alt_2_ra_dec against astropy-computed reference values."""

    def test_absolute_ra_dec_from_az_alt(self):
        """Az/Alt for Altair at OCA 2023-03-27 15:00:16 → recover FK5 J2000 RA/Dec.

        Reference az=304.29692°, alt=37.31958° → RA=297.6945°, Dec=8.8675° (Altair J2000)
        """
        OCA = {'latitude': -24.598056, 'longitude': -70.196389, 'elevation': 2817}
        time = datetime.datetime(2023, 3, 27, 15, 0, 16)
        ra, dec = az_alt_2_ra_dec(
            az=304.29691654663344,
            alt=37.31957620407559,
            latitude=OCA['latitude'],
            longitude=OCA['longitude'],
            elevation=OCA['elevation'],
            time=time,
        )
        self.assertAlmostEqual(ra, 297.6945, places=3)
        self.assertAlmostEqual(dec, 8.8675, places=3)


class TestRadecToJ2000(unittest.TestCase):
    """Test radec_to_j2000 precession accuracy against astropy FK5 reference."""

    def test_j2000_identity(self):
        """Precessing from J2000 to J2000 should return identical coordinates."""
        ra_in, dec_in = 83.82208, -5.39111
        ra_out, dec_out = radec_to_j2000(ra_in, dec_in, epoch='J2000')
        self.assertAlmostEqual(ra_out, ra_in, places=8)
        self.assertAlmostEqual(dec_out, dec_in, places=8)

    def test_j2010_to_j2000(self):
        """Precess Betelgeuse-like coords from FK5 J2010 to J2000.

        astropy reference (FK5 J2010→J2000): ra=83.69919°, dec=-5.39716°
        """
        ra_j2000, dec_j2000 = radec_to_j2000(83.82208, -5.39111, epoch='J2010')
        self.assertAlmostEqual(ra_j2000, 83.69919, places=3)
        self.assertAlmostEqual(dec_j2000, -5.39716, places=3)

    def test_j2015_5_to_j2000(self):
        """Precess Sirius-like coords from FK5 J2015.5 to J2000.

        astropy reference (FK5 J2015.5→J2000): ra=101.11382°, dec=-16.69924°
        """
        ra_j2000, dec_j2000 = radec_to_j2000(101.287, -16.716, epoch='J2015.5')
        self.assertAlmostEqual(ra_j2000, 101.11382, places=3)
        self.assertAlmostEqual(dec_j2000, -16.69924, places=3)

    def test_precession_changes_coords(self):
        """Precessing from J2010 should differ measurably from J2000 input."""
        ra_j2000, _ = radec_to_j2000(83.82208, -5.39111, epoch='J2000')
        ra_j2010, _ = radec_to_j2000(83.82208, -5.39111, epoch='J2010')
        self.assertGreater(abs(ra_j2000 - ra_j2010), 0.01)


class TestAzAlt2RaDecAstropy(unittest.TestCase):

    @unittest.skip("Flaky: uses datetime.now() and requires IERS data download unavailable in CI")
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
