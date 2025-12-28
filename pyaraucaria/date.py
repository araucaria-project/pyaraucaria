"""
Handles conversion to and from julian dates and heliocentric correction

Compatible with python 2.7 and 3.6+

Useful functions
----------------
- `datetime_to_julian`
- `julian_to_datetime`

Licence: MIT
"""

import logging
import re
import datetime
from typing import Union

import numpy as np
from astropy.coordinates import SkyCoord, EarthLocation
from astropy.time import Time
import astropy.units as u
from pyaraucaria.erfa import eraCal2jd

logger = logging.getLogger('dates')

import numpy as __np

_datetime_regexps = [re.compile(r) for r in [
    r'^\s*(?P<Y>\d\d\d\d)[-/\.](?P<M>\d\d)[-/\.](?P<D>\d\d)((?:\s+|T)(?P<h>\d?\d):(?P<m>\d?\d)(?::(?P<s>\d?\d(?:\.\d*)?))?)?(?:Z|(?P<tzs>[+-])(?P<tzh>\d\d):?(?P<tzm>\d\d))?(\s.*)?$',
    r'^\s*(?P<D>\d?\d)[-/](?P<M>\d?\d)[-/](?P<Y>\d{2}(?:\d{2})?)((?:\s+)(?P<h>\d?\d):(?P<m>\d?\d)(?:(?::(?P<s>\d?\d(?:\.\d*)?)))?)?(\s.*)?$',
    r'^\s*(?P<Y>\d\d\d\d)(?P<M>\d\d)(?P<D>\d\d)((?: +|T)?(?P<h>\d\d)(?P<m>\d\d)(?:(?P<s>\d?\d(?:\.\d*)?))?)?(?:Z|(?P<tzs>[+-])(?P<tzh>\d\d):?(?P<tzm>\d\d))?(\s.*)?$',
]]


def julian_to_iso(jd, seconds_decimals=3, separator='T', include_date=True, include_time=True):
    # type: (float, int, str, bool, bool) -> str
    """ Formats ISO 8601-like date time from julian date"""
    duple = julian_to_tuple(jd)
    dstr = tuple_to_iso(*duple, seconds_decimals=seconds_decimals, separator=separator,
                        include_date=include_date, include_time=include_time)
    return dstr

def expand_julian(jd):
    # type: (float) -> float
    """Expands Reduced Julian Day (not the Modified Julian Day!) to full Julian day"""
    if jd < 10000:
        return jd + 2450000
    elif jd < 100000:
        return jd + 2400000
    else:
        return jd

vexpand_julian = __np.vectorize(expand_julian)


def julian_to_tuple(jd):
    # type: (float) -> (int, int, int, int, int, float)

    jd = expand_julian(jd)
    ijd = int(jd + 0.5)
    f = ijd + 1363 + (((4 * ijd + 274277) // 146097) * 3) // 4
    e = 4 * f + 3
    g = (e % 1461) // 4
    h = 5 * g + 2
    D = (h % 153) // 5 + 1
    M = ((h // 153) + 2) % 12 + 1
    Y = e // 1461 - 4716 + (14 - M) // 12
    pd = jd+0.5 - ijd
    ph = 24 * pd
    HH = int(ph)
    pm = (ph % 1) * 60
    MM = int(pm)
    SS = (pm % 1) * 60
    return Y, M, D, HH, MM, SS

def tuple_to_iso(Y, M, D, h, m, s,
                 seconds_decimals=3, separator='T', include_date=True, include_time=True):
    # type: (int, int, int, int, int, float, int, str, bool, bool) -> str

    if include_time:
        seconds_format = "2.0" if seconds_decimals == 0 else str(seconds_decimals+3) + "." + str(seconds_decimals)
        if include_date:
            fstr = "{0:4d}-{1:02d}-{2:02d}{6:}{3:02d}:{4:02d}:{5:0" + seconds_format + "f}"
        else:
            fstr = "{3:02d}:{4:02d}:{5:0" + seconds_format + "f}"
    else:
        fstr = "{0:4d}-{1:02d}-{2:02d}"
    return fstr.format(Y, M, D, h, m, s, separator)


def datetime_to_julian(date):
    # type: (Union[str, datetime.datetime]) -> float
    """ Parse date (with optional time and/or timezone) string to julian date

    If timezone is specified (e.g. +02:00) correction is applied in order to return julian date in UT

    Warning
    -------
    For ambiguous date formats like 'xx/yy/zz' or 'zz-yy-zz' the assumed order is DAY, MONTH, YEAR (DD/MM/YY)!

    Parameters
    ----------
    date : str or datetime.datetime
        Date/time to be parsed, supports most of iso8601 format and some other popular representations

    Returns
    -------
    float
        julian date
    """
    dtuple = datetime_to_tuple(date)
    return tuple_to_julian(*dtuple)

def datetime_to_julian_old(date):
    # type: (Union[str, datetime.datetime]) -> float
    """ Parse date (with optional time and/or timezone) string to julian date

    If timezone is specified (e.g. +02:00) correction is applied in order to return julian date in UT

    Warning
    -------
    For ambiguous date formats like 'xx/yy/zz' or 'zz-yy-zz' the assumed order is DAY, MONTH, YEAR (DD/MM/YY)!

    Parameters
    ----------
    date : str or datetime.datetime
        Date/time to be parsed, supports most of iso8601 format and some other popular representations

    Returns
    -------
    float
        julian date
    """
    dtuple = datetime_to_tuple(date)
    return tuple_to_julian_old(*dtuple)

def tuple_to_julian(Y, M, D, h, m, s, tzs=0, tzh=0, tzm=0):
    # type: (int, int, int, int, int, float, int, int, int) -> (float)
    jd0, jd1, js = eraCal2jd(Y, M, D)
    if js != 0:
        raise ValueError('Date (YMD = {}, {}, {}) not correct'.format(Y, M, D))
    jd3 = h / 24.0 + m / 1440.0 + s / 86400.0 - tzs * tzh / 24.0 - tzs * tzm / 1440.0
    return jd0 + jd1 + jd3

def tuple_to_julian_old(Y, M, D, h, m, s, tzs=0, tzh=0, tzm=0):
    # type: (int, int, int, int, int, float, int, int, int) -> (float)
    MM = int((M - 14) / 12.0)
    JD = (1461 * (Y + 4801 + MM)) // 4 + int((367 * (M - 2 - 12 * (1 + MM))) / 12.0) \
        - (3 * ((Y + 4901 + MM) // 100)) // 4 + D - 32075 + (h - 12) / 24.0 + m / 1440.0 + s / 86400.0\
        - tzs * tzh / 24.0 - tzs * tzm / 1440.0
    return JD

def datetime_to_tuple(date):
    # type: (str) -> (int, int, int, int, int, float, int, int, int)
    """ Parse date (with optional time and/or timezone) string.

    Warning
    -------
    For ambiguous date formats like 'xx/yy/zz' or 'zz-yy-zz' the assumed order is DAY, MONTH, YEAR (DD/MM/YY)!

    Parameters
    ----------
    date : str
        Date/time to be parsed, supports most of iso8601 format and some other popular representations

    Returns
    -------
    (int, int, int, int, int, float, int, int, int)
        Returns tuple (year, month, day, hours, minutes, seconds_with_decimals, tz_sign, tz_hours, tz_minutes)

    Examples
    --------
    >>> datetime_to_tuple('1999-08-12T23:34:1Z')
    (1999, 8, 12, 23, 34, 1.0, 0, 0, 0)
    >>> datetime_to_tuple('1999-08-12T23:34:17.134')
    (1999, 8, 12, 23, 34, 17.134, 0, 0, 0)
    >>> datetime_to_tuple('1999-08-12  23:34:17.134')
    (1999, 8, 12, 23, 34, 17.134, 0, 0, 0)
    >>> datetime_to_tuple('1999-08-12T23:34:17.134Z')
    (1999, 8, 12, 23, 34, 17.134, 0, 0, 0)
    >>> datetime_to_tuple(' 1999-08-12T23:34:17. ')
    (1999, 8, 12, 23, 34, 17.0, 0, 0, 0)
    >>> datetime_to_tuple('1999-08-12T23:34:17')
    (1999, 8, 12, 23, 34, 17.0, 0, 0, 0)
    >>> datetime_to_tuple('1999-08-12T23:34')
    (1999, 8, 12, 23, 34, 0.0, 0, 0, 0)
    >>> datetime_to_tuple('1999-08-12T23:34:17.134-02:30')
    (1999, 8, 12, 23, 34, 17.134, -1, 2, 30)
    >>> datetime_to_tuple('1999-08-12T23:34:17-0230')
    (1999, 8, 12, 23, 34, 17.0, -1, 2, 30)
    >>> datetime_to_tuple('1999-08-12T23:34:17.134-02:21')
    (1999, 8, 12, 23, 34, 17.134, -1, 2, 21)
    >>> datetime_to_tuple('1999-08-12')
    (1999, 8, 12, 0, 0, 0.0, 0, 0, 0)
    >>> datetime_to_tuple('1999-08-12Z')
    (1999, 8, 12, 0, 0, 0.0, 0, 0, 0)
    >>> datetime_to_tuple('1999-08-12-02:21')
    (1999, 8, 12, 0, 0, 0.0, -1, 2, 21)
    >>> datetime_to_tuple('1999/08/12 23:34:17.134')
    (1999, 8, 12, 23, 34, 17.134, 0, 0, 0)
    >>> datetime_to_tuple('05/04/87')
    (1987, 4, 5, 0, 0, 0.0, 0, 0, 0)
    >>> datetime_to_tuple('06/09/01')
    (2001, 9, 6, 0, 0, 0.0, 0, 0, 0)
    >>> datetime_to_tuple('06/09/01 12:33:01.12')
    (2001, 9, 6, 12, 33, 1.12, 0, 0, 0)
    >>> datetime_to_tuple('06/09/01 12:33')
    (2001, 9, 6, 12, 33, 0.0, 0, 0, 0)
    >>> datetime_to_tuple('5/4/87')
    (1987, 4, 5, 0, 0, 0.0, 0, 0, 0)
    >>> datetime_to_tuple('2-12-93')
    (1993, 12, 2, 0, 0, 0.0, 0, 0, 0)
    >>> datetime_to_tuple('19990812')
    (1999, 8, 12, 0, 0, 0.0, 0, 0, 0)
    >>> datetime_to_tuple('19990812 235959')
    (1999, 8, 12, 23, 59, 59.0, 0, 0, 0)
    >>> datetime_to_tuple('19990812T235959')
    (1999, 8, 12, 23, 59, 59.0, 0, 0, 0)
    >>> datetime_to_tuple('19990812235959')
    (1999, 8, 12, 23, 59, 59.0, 0, 0, 0)
    >>> datetime_to_tuple('19990812235959+0400')
    (1999, 8, 12, 23, 59, 59.0, 1, 4, 0)
    >>> datetime_to_tuple('19990812')
    (1999, 8, 12, 0, 0, 0.0, 0, 0, 0)
    """
    if isinstance(date, datetime.datetime):
        return date.year, date.month, date.day, date.hour, date.minute, date.second + date.microsecond / 1000000.0
    for regexp in _datetime_regexps:
        try:
            res = regexp.search(date)
        except TypeError as e:
            logger.error('Rexp date failed for argument of type ', type(date), ' and value: ', date)
            raise e
        if res:
            rg = res.groupdict()
            return tuple(t(rg[v] if rg.get(v, None) is not None else 0) for v, t in [
                ('Y', lambda x: correct_year(int(x)) if x else 0),
                ('M', lambda x: int(x) if x else 0),
                ('D', lambda x: int(x) if x else 0),
                ('h', lambda x: int(x) if x else 0),
                ('m', lambda x: int(x) if x else 0),
                ('s', lambda x: float(x) if x else 0.0),
                ('tzs', lambda x: (-1 if x == '-' else 1) if x else 0),
                ('tzh', lambda x: int(x) if x else 0),
                ('tzm', lambda x: int(x) if x else 0),
            ])
    logger.error('Could no parse this datetime: %s', date)
    raise ValueError('Could no parse this datetime: %s', date)


def helio_corr(jd, ra, dec, longitude=None, latitude=None, elevation=None):
    # type: (float, Union[float, str], Union[float, str], Union[float, None], Union[float, None], Union[float, None]) -> float
    """Calculates heliocentric correction for a given julian day for object at specific coordinates.
    Calculated from the Earth's Geocenter or, if longitude, latitude and elevation are provided, from that location on Earth.

    Parameters
    ----------
    jd : float
        julian day
    ra, dec : float or str
        coordinates of event. If str: "HH:MM:SS", "DD:MM:SS".
        If float: RA is assumed to be in Hours, Dec in Degrees.
    longitude, latitude, elevation : float, optional
        Observer's longitude and latitude in degrees, elevation in meters.
        If not provided, Earth's Geocenter is used - consistent with HJD definition.

    Returns
    -------
    float
        Correction in days (add this to JD to get HJD)
    """

    if longitude is not None and latitude is not None and elevation is not None:
        # Define Observer at specified location
        geocenter = EarthLocation(lat=latitude * u.deg,
                                  lon=longitude * u.deg,
                                  height=elevation * u.m)
    else:
        # 1. Define Observer at Earth's Geocenter (0, 0, 0)
        # This creates a location at the center of the Earth - consistent with HJD definition.
        geocenter = EarthLocation.from_geocentric(0, 0, 0, unit=u.m)

    # 2. Create Time object with location attached
    t = Time(jd, format='jd', scale='utc', location=geocenter)

    # 3. Create Target Coordinate
    # units=(u.hourangle, u.deg) handles both string parsing ("05:05:37")
    # and floats (assumed hours for RA, degrees for Dec) correctly.
    target = SkyCoord(ra, dec, unit=(u.hourangle, u.deg), frame='icrs')

    # 4. Calculate Light Travel Time
    # Astropy calculates vector: Sun -> Earth_Center
    ltt = t.light_travel_time(target, kind='heliocentric')

    return ltt.to_value('day')

def correct_year(year):
    # type: (int) -> int
    return year if year > 100 else (year + 1900 if year > 40 else year + 2000)


def get_oca_jd(jd: Union[float, np.ndarray]) -> Union[float, np.ndarray]:
    """
    Func. calculate oca julian date
    @param jd: Input julian date
    @return: Oca jd date
    """
    return jd - 2460000

