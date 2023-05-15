"""
Handles fast transformation of decimal (float) and sexagesimal (str) coordinates

Compatible with python 2.7 and 3.6+

Useful functions
----------------
- `ra_to_decimal`
- `dec_to_decimal`
- `ra_to_sexagesimal`
- `dec_to_sexagesimal`

Licence: MIT
"""

import re
import ephem
import math


def ra_to_decimal(hms):
    # type: (str) -> float
    return hourangle_to_decimal_deg(hms)


def dec_to_decimal(dms):
    # type: (str) -> float
    return deg_to_decimal_deg(dms)


def ra_to_sexagesimal(deg, sep=':', precision=3):
    # type: (float, str, int) -> str
    return to_hourangle_sexagesimal(deg, sep=sep, sign=False, precision=precision)


def dec_to_sexagesimal(deg, sep=':', precision=3):
    # type: (float, str, int) -> str
    return to_degminsec_sexagesimal(deg, sep=sep, sign=True, precision=precision)


def to_hourangle_sexagesimal(deg, sep=':', sign=False, precision=3):
    # type: (float, str, bool, int) -> str
    return format_sexagesimal(deg, 24.0/360, sign=sign, sep=sep, precision=precision)


def to_degminsec_sexagesimal(deg, sep=':', sign=True, precision=3):
    # type: (float, str, bool, int) -> str
    return format_sexagesimal(deg, 1, sign=sign, sep=sep, precision=precision)


def hourangle_to_decimal_deg(hms):
    # type: (str) -> float
    try:
        val = float(hms)
    except ValueError:
        sign, h, m, s = parse_sexagesimal(hms)
        val = sign * ((((s / 60) + m) / 60) + h) / 24 * 360
    return round(val, 6)


def deg_to_decimal_deg(dms):
    # type: (str) -> float
    try:
        val = float(dms)
    except ValueError:
        sign, d, m, s = parse_sexagesimal(dms)
        val = sign * ((((s / 60) + m) / 60) + d)
    return round(val, 6)

_sexagesimal_parser = re.compile(r"(?P<sign>[+\-])?(?P<A>\d\d?)[ :\-hHₕ°](?P<B>\d\d?)[ :\-mMₘ′'](?P<C>\d\d?(?:\.\d*)?)")


def parse_sexagesimal(sexagesimal):
    # type: (str) -> (int, int, int, float)
    """Returns quadruple: (sign, dec/h, min, sec)"""
    v = _sexagesimal_parser.match(sexagesimal)
    if v is None:
        raise ValueError('{} can not be converted to decimal representation'.format(sexagesimal))
    v = v.groupdict()
    if v['sign'] == '-':
        sign = -1
    else:
        sign = 1
    return sign, int(v['A']), int(v['B']), float(v['C'])


def format_sexagesimal(deg, multiplier, sign, sep=':', precision=3):
    # type: (float, float, bool, str, int) -> str
    if len(sep) == 1:
        sep = sep * 2

    sig = 1 if deg > 0 else -1
    h = sig * multiplier * deg
    m = (h * 60) % 60
    s = (m * 60) % 60
    return ('{}{:02d}{}{:02d}{}{:0'+ f'{precision + 3}'+ '.' + f'{precision}' + 'f}{}').format(
        ('+' if sig > 0 else '-') if sign else '',
        int(h),
        sep[0],
        int(m),
        sep[1],
        s,
        sep[2] if len(sep) == 3 else ''
    )

def ra_dec_epoch(ra: str, dec: str, longitude: str, latitude: str, elevation: str, epoch:str):
    """
    Func calculates coordinates (ra, dec) for site longitude latitude elevation and epoch=now
    Parameters
    ----------
    longitude - site longitude
    latitude - site latitude
    elevation - site elevation
    ra - ra for given epoch
    dec - dec for given epoch
    epoch - epoch of (ra, dec), example: '2000/01/01', '2000'

    Returns
    -------
    (ra, dec) for epoch=now and given site
    """
    site = ephem.Observer()
    site.date = ephem.now()
    site.lon = longitude
    site.lat = latitude
    site.elevation = float(elevation)
    star = ephem.FixedBody()
    star._ra = ephem.hours(ra)
    star._dec = ephem.degrees(dec)
    star._epoch = epoch
    star.compute(site)
    return str(star.g_ra), str(star.g_dec)

def ra_dec_2_az_alt(ra: str, dec: str, longitude: str, latitude: str,
                    elevation: str, epoch:str, time = None):
    """
    Func calculate alt, az from ra, dec for given site and time
    Parameters
    ----------
    ra - ra for given epoch
    dec - dec for given epoch
    longitude - site longitude
    latitude - site latitude
    elevation - site elevation
    epoch - epoch of (ra, dec), example: '2000/01/01', '2000'
    time - time UTC

    Returns
    -------
    (alt, az) for given parameters
    """
    site = ephem.Observer()
    if time:
        site.date = time
    else:
        site.date = ephem.now()
    site.lon = longitude
    site.lat = latitude
    site.elevation = float(elevation)
    star = ephem.FixedBody()
    star._ra = ephem.hours(ra)
    star._dec = ephem.degrees(dec)
    star._epoch = epoch
    star.compute(site)

    return str(star.az), str(star.alt)

def deg_str_to_deg(deg: str):
    """
    Converts degrees from str XX:XX:XX to float
    :param deg: degrees in XX:XX:XX (string) format
    :return: degrees in float
    """
    w = deg.split(":")
    if int(w[0]) < 0:
        sign = -1
    else:
        sign = 1
    return sign * (abs(int(w[0])) + int(w[1])/60 + float(w[2])/3600)

def site_sidereal_time(longitude: str, latitude: str,
                    elevation: str, time = None, deg_output: bool = False) -> str or float:
    """
    Func returns site sidereal time
    Parameters
    ----------
    longitude - site longitude
    latitude - site latitude
    elevation - site elevation
    time - make calculation for different time than now, defalut None
    deg_output - if True degrees output

    Returns
    -------
    site sidereal time (str)
    """
    site = ephem.Observer()
    if time:
        site.date = time
    else:
        site.date = ephem.now()
    site.lon = longitude
    site.lat = latitude
    site.elevation = float(elevation)
    sidereal_time = (site.sidereal_time()).__str__()
    if deg_output:
        sidereal_time = hourangle_to_decimal_deg(sidereal_time.__str__())
    return sidereal_time

def az_alt_2_ra_dec(az: float or str, alt: float or str, longitude: str, latitude: str,
                    elevation: str, time = None, ra_hour: bool = False, dec_sex: bool = False):
    """
    Func calculate alt, az from ra, dec for given site and time
    Parameters
    ----------
    az - azimuth
    alt - altitude
    longitude - site longitude
    latitude - site latitude
    elevation - site elevation
    time - time UTC
    ra_hour - if True returns ra in hourangle sexagesimal, else in degrees float
    dec_sex - if True returns dec in sexagesimal, else in degrees float

    Returns
    -------
    (ra, dec) for given parameters
    """
    site = ephem.Observer()
    if time:
        site.date = time
    else:
        site.date = ephem.now()
    site.lon = longitude
    site.lat = latitude
    site.elevation = float(elevation)

    if isinstance(az, str) and isinstance(alt, str):
        az = az
        alt = alt
    elif isinstance(az, float) and isinstance(alt, float):
        az = math.radians(az)
        alt = math.radians(alt)
    else:
        raise ValueError
    _ra, _dec = site.radec_of(az, alt)
    _ra = math.degrees(_ra)
    _dec = math.degrees(_dec)
    if ra_hour:
        ra = ra_to_sexagesimal(_ra)
    else:
        ra = _ra
    if dec_sex:
        dec = dec_to_sexagesimal(_dec)
    else:
        dec = _dec

    return ra, dec
