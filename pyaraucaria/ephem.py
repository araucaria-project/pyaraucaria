from __future__ import annotations

from datetime import datetime
from enum import Enum
from functools import lru_cache
from typing import Tuple, Sequence

import numpy as np
from astropy.units import u
from astropy.coordinates import EarthLocation, get_moon, AltAz, SkyCoord, get_sun, get_body
from astropy.time import Time

oca_location = EarthLocation.from_geodetic(lon=-70.20128, lat=-24.59867, height=2817)


def get_body_ephem(body: SkyCoord,
                   times: list[datetime] | list[Time] | Time | None = None,
                   location: EarthLocation | Tuple[float, float] | Tuple[float, float, float] | Tuple
                   [float, float, float, str] = oca_location) -> (list[float], list[float], list[float], list[float]):
    """ For given body and list of moments, return lists of ra, dec, az, alt for each moment.
    """
    location = as_location(location)
    times = as_astropy_time(times)
    altaz = AltAz(obstime=times, location=location)
    body_alt_az = body.transform_to(altaz)
    ra = np.atleast_1d(body.ra.deg).tolist()
    dec = np.atleast_1d(body.dec.deg).tolist()
    alt = np.atleast_1d(body_alt_az.alt.deg).tolist()
    az = np.atleast_1d(body_alt_az.az.deg).tolist()
    return ra, dec, az, alt


def get_moon_ephem(times: list[datetime] | list[Time] | Time | None = None,
                   location: EarthLocation | Tuple[float, float] | Tuple[float, float, float] | Tuple
                   [float, float, float, str] = oca_location) -> list[dict[str, float]]:
    """ For given list of moments, return list of dicts with ra, dec, alt, az of the Moon for each moment.

    See Also:
        get_body_ephem
    """
    location = as_location(location)
    times = as_astropy_time(times)
    body = get_moon(times, location)
    return get_body_ephem(body, times, location)


def get_sun_ephem(times: list[datetime] | list[Time] | Time | None = None,
                  location: EarthLocation | Tuple[float, float] | Tuple[float, float, float] | Tuple
                  [float, float, float, str] = oca_location) -> list[dict[str, float]]:
    """ For given list of moments, return list of dicts with ra, dec, alt, az of the Sun for each moment.

    See Also:
        get_body_ephem
    """
    times = as_astropy_time(times)
    body = get_sun(times)
    return get_body_ephem(body, times, location)


class WhichMoments(Enum):
    LAST = 'last'
    NEXT = 'next'
    CLOSEST = 'closest'


def get_body_alt_times(body: SkyCoord, altitudes: Sequence[float],
                       which: str | WhichMoments = WhichMoments.CLOSEST,
                       around_time: datetime | Time | None = None,
                       ) -> list[datetime]:
    """ For given body and list of altitudes, return list of moments when body is at each altitude.

    Args:
        body: SkyCoord of body
        altitudes: list of altitudes in degrees
        which: WhichMoments enum or string 'last', 'next', 'closest'
        around_time: datetime or Time object to search around. If None, uses current time.

    Returns:
        list of datetime objects
    """
    around_time = as_astropy_time(around_time)
    if which == WhichMoments.LAST:
        around_time = around_time - 12 * u.hour
    elif which == WhichMoments.NEXT:
        around_time = around_time + 12 * u.hour


@lru_cache(maxsize=10)
def _get_alt_grid(body: str | SkyCoord, location: EarthLocation,
                 t_start: Time, t_end: Time,
                 ) -> ([float], [float]):
    """ Calculates (and caches) altitudes and azimuth for selected period with minute time resolution"""
    dt = np.arange(0, (t_end - t_start).to(u.minute).value, 1) * u.minute

    times = t_start + dt

    body = as_sky_coord(body, times, location)
    _, _, alt, az = get_body_ephem(body, times, location)




def as_location(location: EarthLocation | str | Tuple[float, float] | Tuple[float, float, float]
                | Tuple[float, float, float, str] = oca_location):
    """Takes a location in various formats and returns an EarthLocation object

    For tuples, the elements are interpreted as (latitude, longitude, elevation=0, ellipsoid='WGS84'),
    str is interpreted as a name of site from astropy.coordinates.EarthLocation.get_site_names()
    """
    if not isinstance(location, EarthLocation):
        if isinstance(location, str):
            location = EarthLocation.of_site(location)
        else:
            labels = ['lon', 'lat', 'height', 'ellipsoid']
            kwargs = dict(zip(labels, location))
            location = EarthLocation.from_geodetic(**kwargs)
    return location

def as_sky_coord(body: str | SkyCoord,
                 times: list[datetime] | list[Time] | Time | None = None,
                 location: EarthLocation | Tuple[float, float] | Tuple[float, float, float] | Tuple = oca_location) \
        -> SkyCoord:
    """Takes a body in various formats and returns a SkyCoord object

    """
    location = as_location(location)
    times = as_astropy_time(times)
    if isinstance(body, str):
        body = get_body(body, times, location)
    elif times is not None:
        body = SkyCoord(body, times, location)
    return body


def as_astropy_time(times: list[datetime] | list[Time] | Time | None = None):
    """Takes a list of datetime objects and returns an astropy.time.Time object

    If times is None, returns current time: `Time.now()`
    """
    if times is None:
        times = Time.now()
    else:
        times = Time(times)
    return times
