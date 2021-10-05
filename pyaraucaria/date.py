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
logger = logging.getLogger('dates')

_datetime_regexps = [re.compile(r) for r in [
    r'(?P<Y>\d\d\d\d)[-/](?P<M>\d\d)[-/](?P<D>\d\d)((?:\s+|T)(?P<h>\d?\d):(?P<m>\d?\d)(?::(?P<s>\d?\d(?:\.\d*)?))?)?(?:Z|(?P<tzs>[+-])(?P<tzh>\d\d):?(?P<tzm>\d\d))?',
    r'(?P<D>\d?\d)[-/](?P<M>\d?\d)[-/](?P<Y>\d?\d)((?:\s+)(?P<h>\d?\d):(?P<m>\d?\d)(?:(?::(?P<s>\d?\d(?:\.\d*)?)))?)?(?P<tzs>)(?P<tzh>)(?P<tzm>)',
    r'(?P<Y>\d\d\d\d)(?P<M>\d\d)(?P<D>\d\d)((?: +|T)?(?P<h>\d\d)(?P<m>\d\d)(?:(?P<s>\d?\d(?:\.\d*)?))?)?(?:Z|(?P<tzs>[+-])(?P<tzh>\d\d):?(?P<tzm>\d\d))?',
]]

def datetime_to_julian(date):
    # type: (str) -> float
    """ Parse date (with optional time and/or timezone) string to julian date

    If timezone is specified (e.g. +02:00) correction is applied in order to return julian date in UT

    Warning
    -------
    For ambiguous date formats like 'xx/yy/zz' or 'zz-yy-zz' the assumed order is DAY, MONTH, YEAR (DD/MM/YY)!

    Parameters
    ----------
    date : str
        Date/time to be parsed, supports most of iso8601 format and some other popular representations

    Returns
    -------
    float
        julian date
    """
    Y, M, D, h, m, s, tzs, tzh, tzm = parse_datetime(date)
    return (1461 * (Y + 4801 + (M - 14) // 12)) // 4 + (367 * (M - 2 - 12 * (1 + (M - 14) // 12))) // 12 \
        - (3 * ((Y + 4901 + (M - 14) // 12) // 100)) // 4 + D - 32075 + (h - 12) / 24.0 + m / 1440.0 + s / 86400.0\
        - tzs * tzh / 24.0 - tzs * tzm / 1440.0


def parse_datetime(date):
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
    >>> parse_datetime('1999-08-12T23:34:1Z')
    (1999, 8, 12, 23, 34, 1.0, 0, 0, 0)
    >>> parse_datetime('1999-08-12T23:34:17.134')
    (1999, 8, 12, 23, 34, 17.134, 0, 0, 0)
    >>> parse_datetime('1999-08-12  23:34:17.134')
    (1999, 8, 12, 23, 34, 17.134, 0, 0, 0)
    >>> parse_datetime('1999-08-12T23:34:17.134Z')
    (1999, 8, 12, 23, 34, 17.134, 0, 0, 0)
    >>> parse_datetime(' 1999-08-12T23:34:17. ')
    (1999, 8, 12, 23, 34, 17.0, 0, 0, 0)
    >>> parse_datetime('1999-08-12T23:34:17')
    (1999, 8, 12, 23, 34, 17.0, 0, 0, 0)
    >>> parse_datetime('1999-08-12T23:34')
    (1999, 8, 12, 23, 34, 0.0, 0, 0, 0)
    >>> parse_datetime('1999-08-12T23:34:17.134-02:30')
    (1999, 8, 12, 23, 34, 17.134, -1, 2, 30)
    >>> parse_datetime('1999-08-12T23:34:17-0230')
    (1999, 8, 12, 23, 34, 17.0, -1, 2, 30)
    >>> parse_datetime('1999-08-12T23:34:17.134-02:21')
    (1999, 8, 12, 23, 34, 17.134, -1, 2, 21)
    >>> parse_datetime('1999-08-12')
    (1999, 8, 12, 0, 0, 0.0, 0, 0, 0)
    >>> parse_datetime('1999-08-12Z')
    (1999, 8, 12, 0, 0, 0.0, 0, 0, 0)
    >>> parse_datetime('1999-08-12-02:21')
    (1999, 8, 12, 0, 0, 0.0, -1, 2, 21)
    >>> parse_datetime('1999/08/12 23:34:17.134')
    (1999, 8, 12, 23, 34, 17.134, 0, 0, 0)
    >>> parse_datetime('05/04/87')
    (1987, 4, 5, 0, 0, 0.0, 0, 0, 0)
    >>> parse_datetime('06/09/01')
    (2001, 9, 6, 0, 0, 0.0, 0, 0, 0)
    >>> parse_datetime('06/09/01 12:33:01.12')
    (2001, 9, 6, 12, 33, 1.12, 0, 0, 0)
    >>> parse_datetime('06/09/01 12:33')
    (2001, 9, 6, 12, 33, 0.0, 0, 0, 0)
    >>> parse_datetime('5/4/87')
    (1987, 4, 5, 0, 0, 0.0, 0, 0, 0)
    >>> parse_datetime('2-12-93')
    (1993, 12, 2, 0, 0, 0.0, 0, 0, 0)
    >>> parse_datetime('19990812')
    (1999, 8, 12, 0, 0, 0.0, 0, 0, 0)
    >>> parse_datetime('19990812 235959')
    (1999, 8, 12, 23, 59, 59.0, 0, 0, 0)
    >>> parse_datetime('19990812T235959')
    (1999, 8, 12, 23, 59, 59.0, 0, 0, 0)
    >>> parse_datetime('19990812235959')
    (1999, 8, 12, 23, 59, 59.0, 0, 0, 0)
    >>> parse_datetime('19990812235959+0400')
    (1999, 8, 12, 23, 59, 59.0, 1, 4, 0)
    >>> parse_datetime('19990812')
    (1999, 8, 12, 0, 0, 0.0, 0, 0, 0)
    """
    for regexp in _datetime_regexps:
        res = regexp.search(date)
        if res:
            rg = res.groupdict()
            return tuple(t(rg[v] if rg[v] is not None else 0) for v, t in [
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


def correct_year(year):
    # type: (int) -> int
    return year if year > 100 else (year + 1900 if year > 40 else year + 2000)
