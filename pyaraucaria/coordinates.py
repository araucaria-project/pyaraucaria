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
import datetime
from astropy.coordinates import SkyCoord, EarthLocation, FK5
import astropy.units as u


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

_sexagesimal_parser = re.compile(r"(?P<sign>[+\-])?(?P<A>\d\d?\d?)[ :\-hHₕ°](?P<B>\d\d?)[ :\-mMₘ′'](?P<C>\d\d?(?:\.\d*)?)")


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


def ra_dec_epoch(ra, dec, epoch):
    # type: (float, float, str) -> (float, float)
    """
    Func calculates coordinates (ra, dec) for site longitude latitude elevation and epoch
    Parameters
    ----------
    ra - ra for given epoch
    dec - dec for given epoch
    epoch - epoch of param: ra, dec; example: '2000/01/01', '2000'

    Returns
    -------
    (ra, dec) for epoch now
    """
    ob = ephem.Observer()
    ob.epoch = ephem.now()
    star = ephem.FixedBody()
    star._epoch = epoch
    star._ra = math.radians(ra)
    star._dec = math.radians(dec)
    star.compute(ob)

    return math.degrees(star.g_ra), math.degrees(star.g_dec)


def ra_dec_2_az_alt(ra, dec, longitude, latitude, elevation, epoch, time=None):
    # type: (float, float, float, float, float, str, datetime or None) -> (float, float)
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
    site.lon = math.radians(longitude)
    site.lat = math.radians(latitude)
    site.elevation = float(elevation)
    star = ephem.FixedBody()
    star._ra = math.radians(ra)
    star._dec = math.radians(dec)
    star._epoch = epoch
    star.compute(site)

    return math.degrees(star.az), math.degrees(star.alt)


def site_sidereal_time(longitude, latitude, elevation, time=None):
    # type: (float, float, float, datetime or None) -> float
    """
    Func returns site sidereal time
    Parameters
    ----------
    longitude - site longitude
    latitude - site latitude
    elevation - site elevation
    time - make calculation for different time than now, defalut None

    Returns
    -------
    site sidereal time (str)
    """
    site = ephem.Observer()
    if time:
        site.date = time
    else:
        site.date = ephem.now()
    site.lon = math.radians(longitude)
    site.lat = math.radians(latitude)
    site.elevation = elevation

    return hourangle_to_decimal_deg(site.sidereal_time().__str__())


def az_alt_2_ra_dec(az, alt, longitude, latitude, elevation, time = None):
    # type: (float, float, float, float, float, datetime or None) -> (float, float)
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

    Returns
    -------
    (ra, dec) for given parameters
    """
    site = ephem.Observer()
    if time:
        site.date = time
    else:
        site.date = ephem.now()
    site.lon = math.radians(longitude)
    site.lat = math.radians(latitude)
    site.elevation = elevation
    az = math.radians(az)
    alt = math.radians(alt)
    _ra, _dec = site.radec_of(az, alt)

    return math.degrees(_ra), math.degrees(_dec)

def az_alt_2_ra_dec_astropy(az, alt, longitude, latitude, elevation, epoch, calc_time=None):
    # type: (float, float, float, float, float, str, datetime or None) -> (float, float)
    """
    Func calculate alt, az from ra, dec for given site and time
    Parameters
    ----------
    az - azimuth
    alt - altitude
    longitude - site longitude
    latitude - site latitude
    elevation - site elevation
    epoch - ra, dec epoch calculation. Exemple: 'J2000'
    time - time UTC

    Returns
    -------
    (ra, dec) for given parameters
    """
    if calc_time:
        time_s = calc_time
    else:
        time_s = datetime.datetime.now()
    loc = EarthLocation.from_geodetic(lon=longitude * u.degree, lat=latitude * u.degree, height=elevation)
    alt_az = SkyCoord(alt=alt * u.degree, az=az * u.degree, frame='altaz', obstime=time_s, location=loc)
    radec = alt_az.transform_to('icrs')
    radec_equinox = radec.transform_to(FK5(equinox=epoch))

    return radec_equinox.ra.deg, radec_equinox.dec.deg
