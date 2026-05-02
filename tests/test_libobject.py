"""
Tests for date/time utilities in pyaraucaria.libobject.

Reference Julian dates used throughout:
  J2000.0  = 2000-01-01T12:00:00 UTC  → JD 2451545.0  (exact, IAU definition)
  2023-03-27 00:00:00 UTC              → JD 2460030.5
  2023-03-27 15:00:16 UTC              → JD 2460031.125185...
  1980-05-30 18:00:00 UTC (= 1980/05/30.75 in ephem notation) → JD 2444390.25
"""
import datetime
import unittest

from pyaraucaria.libobject import (
    _parse_date_to_time,
    get_jday_now,
    get_jday_date,
    get_ut_now,
    get_ut_date,
    get_lt_date,
    ut2jday,
    utvec2jday,
)


class TestParseDateToTime(unittest.TestCase):
    """Test _parse_date_to_time parses various date formats correctly."""

    def test_naive_datetime_treated_as_utc(self):
        """A naive datetime should be treated as UTC and give the correct JD."""
        dt = datetime.datetime(2000, 1, 1, 12, 0, 0)
        t = _parse_date_to_time(dt)
        self.assertAlmostEqual(t.jd, 2451545.0, places=6)

    def test_aware_datetime_preserved(self):
        """A timezone-aware datetime should round-trip to the same UTC JD."""
        dt = datetime.datetime(2000, 1, 1, 12, 0, 0, tzinfo=datetime.timezone.utc)
        t = _parse_date_to_time(dt)
        self.assertAlmostEqual(t.jd, 2451545.0, places=6)

    def test_ephem_style_datetime_string(self):
        """Standard ephem date string 'YYYY/MM/DD HH:MM:SS' → correct JD."""
        t = _parse_date_to_time('2000/01/01 12:00:00')
        self.assertAlmostEqual(t.jd, 2451545.0, places=6)

    def test_ephem_style_date_only(self):
        """Date-only string 'YYYY/MM/DD' → start of day (00:00:00 UTC)."""
        t = _parse_date_to_time('2023/03/27')
        self.assertAlmostEqual(t.jd, 2460030.5, places=6)

    def test_fractional_day_notation(self):
        """Fractional day 'YYYY/MM/DD.frac' → correct JD.

        1980/05/30.75 = 1980-05-30 18:00:00 UTC → JD 2444390.25
        """
        t = _parse_date_to_time('1980/05/30.75')
        self.assertAlmostEqual(t.jd, 2444390.25, places=5)

    def test_fractional_day_zero(self):
        """Fractional day with .0 suffix → midnight UTC."""
        t = _parse_date_to_time('2023/03/27.0')
        self.assertAlmostEqual(t.jd, 2460030.5, places=6)


class TestGetJdayDate(unittest.TestCase):
    """Test get_jday_date returns correct Julian dates for known reference epochs."""

    def test_j2000_epoch(self):
        """J2000.0 epoch: '2000/01/01 12:00:00' → JD=2451545.0 (exact by definition)."""
        jd = get_jday_date('2000/01/01 12:00:00')
        self.assertAlmostEqual(jd, 2451545.0, places=6)

    def test_known_date_midnight(self):
        """'2023/03/27 00:00:00' → JD=2460030.5."""
        jd = get_jday_date('2023/03/27 00:00:00')
        self.assertAlmostEqual(jd, 2460030.5, places=6)

    def test_fractional_day(self):
        """'1980/05/30.75' (ephem fractional-day notation) → JD=2444390.25."""
        jd = get_jday_date('1980/05/30.75')
        self.assertAlmostEqual(jd, 2444390.25, places=5)

    def test_date_only(self):
        """'2023/03/27' (date-only string) → start-of-day JD=2460030.5."""
        jd = get_jday_date('2023/03/27')
        self.assertAlmostEqual(jd, 2460030.5, places=6)


class TestGetJdayNow(unittest.TestCase):
    """Test get_jday_now returns a plausible current Julian date."""

    def test_jday_now_is_reasonable(self):
        """Current JD should be in the range 2020–2100."""
        jd = get_jday_now()
        self.assertGreater(jd, 2458849.5)   # after 2020-01-01
        self.assertLess(jd, 2816788.5)      # before 2100-01-01

    def test_jday_now_is_float(self):
        """Return type must be a float (or float-compatible)."""
        jd = get_jday_now()
        self.assertIsInstance(float(jd), float)


class TestUt2Jday(unittest.TestCase):
    """Test ut2jday with known reference values."""

    def test_j2000_datetime(self):
        """datetime(2000,1,1,12) → JD=2451545.0 (J2000.0)."""
        dt = datetime.datetime(2000, 1, 1, 12, 0, 0)
        self.assertAlmostEqual(ut2jday(dt), 2451545.0, places=6)

    def test_known_date_time(self):
        """datetime(2023,3,27,15,0,16) → JD≈2460031.1252."""
        dt = datetime.datetime(2023, 3, 27, 15, 0, 16)
        self.assertAlmostEqual(ut2jday(dt), 2460031.125185, places=4)

    def test_midnight_gives_half_jd(self):
        """A midnight datetime gives a JD ending in .5 (Julian day starts at noon)."""
        dt = datetime.datetime(2023, 3, 27, 0, 0, 0)
        jd = ut2jday(dt)
        self.assertAlmostEqual(jd % 1, 0.5, places=6)


class TestUtvec2Jday(unittest.TestCase):
    """Test utvec2jday converts a sequence of datetimes to an array of JDs."""

    def test_two_known_dates(self):
        """Vector of two datetimes → array of two correct JDs."""
        dts = [
            datetime.datetime(2000, 1, 1, 12, 0, 0),
            datetime.datetime(2023, 3, 27, 0, 0, 0),
        ]
        jds = utvec2jday(dts)
        self.assertAlmostEqual(float(jds[0]), 2451545.0, places=6)
        self.assertAlmostEqual(float(jds[1]), 2460030.5, places=6)

    def test_length_preserved(self):
        """Output array has the same length as the input sequence."""
        dts = [datetime.datetime(2000, 1, 1, 12) for _ in range(5)]
        jds = utvec2jday(dts)
        self.assertEqual(len(jds), 5)


class TestGetUtDate(unittest.TestCase):
    """Test get_ut_date parses ephem-style strings to UTC datetime objects."""

    def test_j2000_string(self):
        """'2000/01/01 12:00:00' → UTC datetime(2000,1,1,12,0,0)."""
        ut = get_ut_date('2000/01/01 12:00:00')
        self.assertEqual(ut.year, 2000)
        self.assertEqual(ut.month, 1)
        self.assertEqual(ut.day, 1)
        self.assertEqual(ut.hour, 12)
        self.assertEqual(ut.minute, 0)
        self.assertEqual(ut.tzinfo, datetime.timezone.utc)

    def test_known_date_time_string(self):
        """'2023/03/27 15:00:16' → correct UTC datetime fields."""
        ut = get_ut_date('2023/03/27 15:00:16')
        self.assertEqual(ut.year, 2023)
        self.assertEqual(ut.month, 3)
        self.assertEqual(ut.day, 27)
        self.assertEqual(ut.hour, 15)
        self.assertEqual(ut.second, 16)
        self.assertEqual(ut.tzinfo, datetime.timezone.utc)

    def test_returns_utc_aware_datetime(self):
        """Returned datetime should always be UTC-aware."""
        ut = get_ut_date('2023/03/27 00:00:00')
        self.assertIsNotNone(ut.tzinfo)
        self.assertEqual(ut.utcoffset(), datetime.timedelta(0))


class TestGetUtNow(unittest.TestCase):
    """Test get_ut_now returns a UTC-aware datetime near the current time."""

    def test_returns_utc_aware(self):
        """Returned datetime should be timezone-aware and UTC."""
        now = get_ut_now()
        self.assertIsNotNone(now.tzinfo)
        self.assertEqual(now.utcoffset(), datetime.timedelta(0))

    def test_year_is_reasonable(self):
        """Year should be in the plausible range for tests running near 2026."""
        now = get_ut_now()
        self.assertGreater(now.year, 2020)
        self.assertLess(now.year, 2100)


class TestGetLtDate(unittest.TestCase):
    """Test get_lt_date applies timezone offsets correctly to a UTC date string."""

    def test_negative_offset(self):
        """UTC midnight shifted by -4 h → previous day 20:00 in UTC-4."""
        lt = get_lt_date('2023/03/27 00:00:00', force_tz=-4)
        # 2023-03-27 00:00 UTC  −  4 h  =  2023-03-26 20:00 (UTC-4)
        self.assertEqual(lt.year, 2023)
        self.assertEqual(lt.month, 3)
        self.assertEqual(lt.day, 26)
        self.assertEqual(lt.hour, 20)
        self.assertEqual(lt.minute, 0)

    def test_positive_fractional_offset(self):
        """UTC 2023-03-27 00:00 shifted by +5.5 h → 2023-03-27 05:30 (IST)."""
        lt = get_lt_date('2023/03/27 00:00:00', force_tz=5.5)
        self.assertEqual(lt.year, 2023)
        self.assertEqual(lt.month, 3)
        self.assertEqual(lt.day, 27)
        self.assertEqual(lt.hour, 5)
        self.assertEqual(lt.minute, 30)

    def test_zero_offset(self):
        """force_tz=0 → same as UTC, date/time fields unchanged."""
        lt = get_lt_date('2023/03/27 12:00:00', force_tz=0)
        self.assertEqual(lt.year, 2023)
        self.assertEqual(lt.month, 3)
        self.assertEqual(lt.day, 27)
        self.assertEqual(lt.hour, 12)


if __name__ == '__main__':
    unittest.main()
