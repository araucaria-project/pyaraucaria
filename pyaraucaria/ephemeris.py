from typing import List, Dict, Union, Optional
import sys
import csv
import warnings
from datetime import datetime, timezone
from functools import lru_cache
import argparse

import numpy as np
from scipy.interpolate import UnivariateSpline

from astropy import units as u
from astropy.time import Time
from astropy.coordinates import SkyCoord, EarthLocation, AltAz, get_body, FK5, get_sun
from astropy.coordinates.errors import UnknownSiteException

from astroplan import moon_illumination
from astroplan import Observer

from pyaraucaria.coordinates import ra_to_decimal, dec_to_decimal


# Local helper: normalize various time-like inputs to an astropy Time
def _ensure_time(t: Optional[Union[Time, datetime]], wrap_scalar: bool = False, **time_kwargs) -> Time:
    """
    Normalize various time-like inputs to an astropy ``Time``.

    Parameters
    ----------
    t : None | Time | datetime | array-like
        Input time value(s).
    wrap_scalar : bool
        If True and the resulting ``Time`` is scalar, return a
        1-element ``Time`` array instead (useful when callers expect
        indexable time arrays). Default False.
    **time_kwargs
        Passed directly to ``astropy.time.Time(...)`` when constructing
        a Time from non-Time inputs (e.g. ``format='jd'``, ``location=...``).
    """
    if t is None:
        out = Time.now()
    elif isinstance(t, Time):
        out = t
    else:
        out = Time(t, **time_kwargs)

    if wrap_scalar and getattr(out, 'isscalar', False):
        # Re-create as a 1-element Time array using the original input
        # (preserving any kwargs such as format/scale/location).
        if isinstance(t, Time):
            return Time([t])
        return Time([t], **time_kwargs)

    return out




# ==========================================
# Module-level cache for next-rise/set queries
# ==========================================

# Tuned for live-poll workloads (TOI etc.) where the same site is queried
# repeatedly at slightly different ``now`` values. With 60-second time
# rounding a single site/altitude/direction occupies one entry per minute.
_NEXT_EVENT_CACHE_SIZE = 256
_TIME_GRAN_SEC = 60.0


def _round_mjd_to_seconds(mjd: float, gran_sec: float = _TIME_GRAN_SEC) -> float:
    """Round an MJD to the nearest ``gran_sec`` boundary (default: 60 s)."""
    sec_per_day = 86400.0
    return round(mjd * sec_per_day / gran_sec) * gran_sec / sec_per_day


@lru_cache(maxsize=_NEXT_EVENT_CACHE_SIZE)
def _next_event_cached(body_kind: str, lat_deg: float, lon_deg: float,
                       height_m: float, altitude_deg: float, direction: str,
                       start_mjd: float) -> Optional[datetime]:
    """LRU-cached next-rise/set computation.

    Keyed on already-rounded scalar inputs so equivalent calls hit the
    cache without floating-point identity issues. Hot path is the
    body-specific :meth:`CelestialBody._alt_grid_deg` computation
    underneath :meth:`CelestialBody._next_event_uncached`.
    """
    location = EarthLocation(
        lat=lat_deg * u.deg, lon=lon_deg * u.deg, height=height_m * u.m)
    start = Time(start_mjd, format='mjd', scale='utc')
    if body_kind == 'sun':
        body: 'CelestialBody' = Sun(location)
    elif body_kind == 'moon':
        body = Moon(location)
    else:
        raise ValueError(f"unknown body_kind {body_kind!r}")
    return body._next_event_uncached(altitude_deg, direction, start)


def clear_next_event_cache() -> None:
    """Drop all cached next-rise/set results (mostly useful in tests)."""
    _next_event_cached.cache_clear()


# ==========================================
# Core Library Classes
# ==========================================

class CelestialBody:
    """
    Base class for celestial calculations handling location and time grids.
    """

    def __init__(self, location: EarthLocation):
        """
        Initialize with an observer's location.

        :param location: astropy.coordinates.EarthLocation object
        """
        self.location = location

    def _get_time_grid(self, start_time: Time, duration_hours: int = 24, step_minutes: int = 5) -> Time:
        """
        Generates a dense time grid for finding altitude crossings using interpolation.
        """
        times = start_time + np.arange(0, duration_hours * 60, step_minutes) * u.min
        return times

    def _find_altitude_crossings(self, target_alt_deg: float, time_grid: Time, alt_grid: u.Quantity) -> Time:
        """
        Numerically finds exact times when the object crosses the target altitude
        using cubic spline interpolation.

        :param target_alt_deg: Target altitude in degrees
        :param time_grid: Grid of astropy Times
        :param alt_grid: Grid of altitudes corresponding to time_grid
        :return: astropy.time.Time object containing crossing moments
        """
        # Shift function so that target altitude becomes zero (root finding)
        alt_values = alt_grid.deg - target_alt_deg
        time_mjd = time_grid.mjd

        # Create spline (s=0 forces the spline to pass through all points)
        spline = UnivariateSpline(time_mjd, alt_values, s=0)

        # Find roots
        roots_mjd = spline.roots()

        # Filter roots to ensure they are within the checked time range
        valid_roots = roots_mjd[(roots_mjd >= time_mjd[0]) & (roots_mjd <= time_mjd[-1])]

        return Time(valid_roots, format='mjd', scale='utc')

    def get_events_by_altitude(self, altitudes_deg: List[float], start_time: Optional[Union[Time, datetime]] = None):
        """Abstract method to be implemented by children."""
        raise NotImplementedError

    def get_ephemeris(self, times: Union[List[Time], List[datetime], Time, datetime]):
        """Abstract method to be implemented by children."""
        raise NotImplementedError

    # Subclasses set this to a short string ('sun', 'moon') to opt into
    # the module-level result cache in :func:`get_next_event_by_altitude`.
    # Bodies whose state is not fully captured by ``self.location``
    # (e.g. ``Star`` with per-instance RA/Dec) leave this ``None`` and
    # fall through to the uncached path.
    _BODY_KIND: Optional[str] = None

    # Step (minutes) of the coarse altitude grid used by the rise/set
    # primitive. The cubic-spline root finder over a 24 h window is
    # accurate to <1 s at horizon crossings even at 10-min spacing,
    # because altitude is near-linear in time at the horizon.
    _NEXT_EVENT_STEP_MIN = 10

    def _alt_grid_deg(self, times: Time) -> np.ndarray:
        """Altitude (degrees) of the body at each grid time.

        Subclass override is the sole hot path for
        :meth:`get_next_event_by_altitude` — keep it lean (one frame
        transform, no extra ICRS/exact-time recomputation).
        """
        raise NotImplementedError

    def _next_event_uncached(self, altitude_deg: float, direction: str,
                             start_time: Time) -> Optional[datetime]:
        times_grid = self._get_time_grid(
            start_time, step_minutes=self._NEXT_EVENT_STEP_MIN)
        alt_deg = np.asarray(self._alt_grid_deg(times_grid), dtype=float)

        time_mjd = times_grid.mjd
        spline = UnivariateSpline(time_mjd, alt_deg - altitude_deg, s=0)
        roots = spline.roots()
        roots = roots[(roots >= time_mjd[0]) & (roots <= time_mjd[-1])]
        if len(roots) == 0:
            return None

        # Direction is inferred from the alternation pattern starting at
        # the first grid point — no extra ephemeris call needed since
        # alt_deg[0] is alt at start_time itself.
        state_above = bool(alt_deg[0] > altitude_deg)
        want_rising = (direction == 'rising')
        for r in roots:
            is_rising = not state_above
            if is_rising == want_rising:
                return Time(float(r), format='mjd',
                            scale='utc').to_datetime(timezone.utc)
            state_above = not state_above
        return None

    def get_next_event_by_altitude(self, altitude_deg: float, direction: str,
                                   start_time: Optional[Union[Time, datetime]] = None
                                   ) -> Optional[datetime]:
        """First crossing of ``altitude_deg`` after ``start_time`` going
        in the requested vertical direction.

        Generic primitive used by :py:func:`calculate_sun_rise_set` and
        :py:func:`calculate_moon_rise_set`; any caller that wants e.g.
        "next moonset at horizon" should use this rather than
        re-implementing the rise/set selection.

        For ``Sun`` and ``Moon`` results are LRU-cached at the module
        level (see :data:`_NEXT_EVENT_CACHE_SIZE`). Inputs are rounded
        to a stable key — location to ~1 µdeg, altitude to 1 µdeg,
        ``start_time`` to a 60-second grid — so repeated calls within
        a minute (e.g. live-poll workloads) hit the cache after the
        first computation. The 60-second rounding shifts the *start*
        of the search window by at most 30 s, which is irrelevant for
        a 24 h next-rise/set query.

        :param altitude_deg: target altitude in degrees
        :param direction: ``'rising'`` (object going up through the
            altitude) or ``'setting'`` (going down).
        :param start_time: when to start searching; defaults to now.
        :return: UTC ``datetime`` of the first matching crossing within
            the next 24 h, or ``None`` if no such crossing exists in
            that window (e.g. polar geometry, or an altitude the body
            never reaches at this site/season). The return type is
            always either ``datetime`` or ``None`` — never an
            astropy masked array.
        """
        if direction not in ('rising', 'setting'):
            raise ValueError(
                f"direction must be 'rising' or 'setting', got {direction!r}")

        tm = _ensure_time(start_time)
        if self._BODY_KIND is None:
            return self._next_event_uncached(float(altitude_deg), direction, tm)

        # Cache key: rounded location + altitude + direction + rounded start.
        loc = self.location
        lat = round(float(loc.lat.deg), 6)
        lon = round(float(loc.lon.deg), 6)
        height = round(float(loc.height.to_value(u.m)), 3)
        alt = round(float(altitude_deg), 6)
        start_mjd = _round_mjd_to_seconds(float(tm.utc.mjd))
        return _next_event_cached(self._BODY_KIND, lat, lon, height,
                                  alt, direction, start_mjd)


class Sun(CelestialBody):
    """
    Class for Solar calculations.
    """

    _BODY_KIND = 'sun'

    def _alt_grid_deg(self, times: Time) -> np.ndarray:
        # ``get_sun`` is the analytic low-precision (arcmin) Sun
        # position — orders of magnitude cheaper than ``get_body('sun',
        # location=...)`` because it skips JPL/built-in ephemerides.
        # The position error translates to ~few-second timing error at
        # rise/set, far below the unmodeled refraction (~35 arcmin at
        # the horizon) and well within what users care about for
        # next-rise/set queries.
        frame = AltAz(obstime=times, location=self.location)
        return get_sun(times).transform_to(frame).alt.deg

    def get_events_by_altitude(self, altitudes_deg: List[float], start_time: Optional[Union[Time, datetime]] = None) -> \
    List[Dict]:
        tm = _ensure_time(start_time)

        # Generate coarse grid for the next 24 hours
        times_grid = self._get_time_grid(tm)
        frame = AltAz(obstime=times_grid, location=self.location)

        # Changed get_sun to get_body('sun', ...) for consistency/compatibility
        sun_coo = get_body("sun", times_grid, location=self.location).transform_to(frame)

        results = []

        for alt_target in altitudes_deg:
            # Find precise times
            event_times = self._find_altitude_crossings(alt_target, times_grid, sun_coo.alt)

            if len(event_times) > 0:
                # Calculate precise coordinates at the found times
                frame_exact = AltAz(obstime=event_times, location=self.location)
                exact_coo_altaz = get_body("sun", event_times, location=self.location).transform_to(frame_exact)
                exact_coo_icrs = get_body("sun", event_times, location=self.location)

                for i, t in enumerate(event_times):
                    results.append({
                        'target_alt': alt_target,
                        # FORCE UTC DATETIME
                        'time_utc': t.to_datetime(timezone.utc),
                        'actual_alt': exact_coo_altaz.alt.deg[i],
                        'az': exact_coo_altaz.az.deg[i],
                        'ra': exact_coo_icrs.ra.deg[i],
                        'dec': exact_coo_icrs.dec.deg[i],
                        'body': 'Sun'
                    })

        results.sort(key=lambda x: x['time_utc'])
        return results

    def get_ephemeris(self, times: Union[List[Time], List[datetime], Time, datetime]) -> List[Dict]:
        t = _ensure_time(times, wrap_scalar=True)

        frame = AltAz(obstime=t, location=self.location)
        sun_altaz = get_body("sun", t, location=self.location).transform_to(frame)
        sun_icrs = get_body("sun", t, location=self.location)

        data = []
        for i in range(len(t)):
            data.append({
                # FORCE UTC DATETIME
                'time_utc': t[i].to_datetime(timezone.utc),
                'alt': sun_altaz[i].alt.deg,
                'az': sun_altaz[i].az.deg,
                'ra': sun_icrs[i].ra.deg,
                'dec': sun_icrs[i].dec.deg,
                'body': 'Sun'
            })
        return data


class Moon(CelestialBody):
    """
    Class for Lunar calculations, including phase.
    """

    _BODY_KIND = 'moon'

    def _alt_grid_deg(self, times: Time) -> np.ndarray:
        # Moon parallax is up to ~1° at the horizon, so we cannot use a
        # geocentric shortcut here — the topocentric ``get_body('moon',
        # location=...)`` path is required for correct rise/set timing.
        frame = AltAz(obstime=times, location=self.location)
        return get_body('moon', times, location=self.location).transform_to(frame).alt.deg

    def get_events_by_altitude(self, altitudes_deg: List[float], start_time: Optional[Union[Time, datetime]] = None) -> \
    List[Dict]:
        tm = _ensure_time(start_time)

        times_grid = self._get_time_grid(tm)
        frame = AltAz(obstime=times_grid, location=self.location)

        # Changed get_moon to get_body('moon', ...)
        moon_coo = get_body("moon", times_grid, location=self.location).transform_to(frame)

        results = []
        for alt_target in altitudes_deg:
            event_times = self._find_altitude_crossings(alt_target, times_grid, moon_coo.alt)

            if len(event_times) > 0:
                frame_exact = AltAz(obstime=event_times, location=self.location)
                exact_coo_altaz = get_body("moon", event_times, location=self.location).transform_to(frame_exact)
                exact_coo_icrs = get_body("moon", event_times, location=self.location)
                phases = moon_illumination(event_times)

                for i, t in enumerate(event_times):
                    results.append({
                        'target_alt': alt_target,
                        # FORCE UTC DATETIME
                        'time_utc': t.to_datetime(timezone.utc),
                        'actual_alt': exact_coo_altaz.alt.deg[i],
                        'az': exact_coo_altaz.az.deg[i],
                        'ra': exact_coo_icrs.ra.deg[i],
                        'dec': exact_coo_icrs.dec.deg[i],
                        'phase': phases[i],
                        'body': 'Moon'
                    })
        results.sort(key=lambda x: x['time_utc'])
        return results

    def get_ephemeris(self, times: Union[List[Time], List[datetime], Time, datetime]) -> List[Dict]:
        t = _ensure_time(times, wrap_scalar=True)

        frame = AltAz(obstime=t, location=self.location)
        moon_altaz = get_body("moon", t, location=self.location).transform_to(frame)
        moon_icrs = get_body("moon", t, location=self.location)
        phases = moon_illumination(t)

        data = []
        for i in range(len(t)):
            data.append({
                # FORCE UTC DATETIME
                'time_utc': t[i].to_datetime(timezone.utc),
                'alt': moon_altaz[i].alt.deg,
                'az': moon_altaz[i].az.deg,
                'ra': moon_icrs[i].ra.deg,
                'dec': moon_icrs[i].dec.deg,
                'phase': phases[i],
                'body': 'Moon'
            })
        return data

    def get_phase(self, times: Union[List[Time], List[datetime], Time, datetime]) -> List[Dict]:
        """
        Optimized method to calculate ONLY Moon phase (illumination fraction).
        Significantly faster than get_ephemeris as it skips coordinate transformations.

        :return: List of dicts [{'time_utc': ..., 'phase': 0.0-1.0, 'body': 'Moon'}]
        """
        # 1. Standardize Input (Same logic as other methods)
        t = _ensure_time(times, wrap_scalar=True)

        # 2. Fast Calculation
        # moon_illumination uses Geocentric coordinates, skipping the heavy
        # Topocentric (AltAz) transformation matrix.
        phases = moon_illumination(t)

        # 3. Format Output
        data = []
        for i in range(len(t)):
            data.append({
                'time_utc': t[i].to_datetime(timezone.utc),
                'phase': phases[i],
                'body': 'Moon'
            })
        return data


class Star(CelestialBody):
    """
    Class for a single Star. RA/Dec are constant.
    """

    def __init__(self, location: EarthLocation, ra: Union[float, str], dec: Union[float, str], name: str = "Star"):
        super().__init__(location)
        ra_deg = ra_to_decimal(ra)
        dec_deg = dec_to_decimal(dec)
        self.coord = SkyCoord(ra=ra_deg * u.deg, dec=dec_deg * u.deg, frame='icrs')
        self.name = name

    def _alt_grid_deg(self, times: Time) -> np.ndarray:
        frame = AltAz(obstime=times, location=self.location)
        return self.coord.transform_to(frame).alt.deg

    def get_events_by_altitude(self, altitudes_deg: List[float], start_time: Optional[Union[Time, datetime]] = None) -> \
    List[Dict]:
        start_time = _ensure_time(start_time)

        times_grid = self._get_time_grid(start_time)
        frame = AltAz(obstime=times_grid, location=self.location)
        # Transform constant SkyCoord to AltAz over the time grid
        star_altaz = self.coord.transform_to(frame)

        results = []
        for alt_target in altitudes_deg:
            event_times = self._find_altitude_crossings(alt_target, times_grid, star_altaz.alt)

            if len(event_times) > 0:
                frame_exact = AltAz(obstime=event_times, location=self.location)
                exact_coo = self.coord.transform_to(frame_exact)

                for i, t in enumerate(event_times):
                    results.append({
                        'target_alt': alt_target,
                        # FORCE UTC DATETIME
                        'time_utc': t.to_datetime(timezone.utc),
                        'actual_alt': exact_coo.alt.deg[i],
                        'az': exact_coo.az.deg[i],
                        'ra': self.coord.ra.deg,
                        'dec': self.coord.dec.deg,
                        'body': self.name
                    })
        results.sort(key=lambda x: x['time_utc'])
        return results

    def get_ephemeris(self, times: Union[List[Time], List[datetime], Time, datetime]) -> List[Dict]:
        t = _ensure_time(times, wrap_scalar=True)

        frame = AltAz(obstime=t, location=self.location)
        star_altaz = self.coord.transform_to(frame)

        data = []
        for i in range(len(t)):
            data.append({
                # FORCE UTC DATETIME
                'time_utc': t[i].to_datetime(timezone.utc),
                'alt': star_altaz[i].alt.deg,
                'az': star_altaz[i].az.deg,
                'ra': self.coord.ra.deg,
                'dec': self.coord.dec.deg,
                'body': self.name
            })
        return data


class Stars(CelestialBody):
    """
    Class for multiple stars, using vectorized calculations.
    """

    def __init__(self, location: EarthLocation, stars_list: List[Dict]):
        """
        :param stars_list: List of dicts [{'id': 'Sirius', 'ra': 101.28, 'dec': -16.7}, ...]
        """
        super().__init__(location)
        self.ids = [s['id'] for s in stars_list]
        ras = [s['ra'] for s in stars_list] * u.deg
        decs = [s['dec'] for s in stars_list] * u.deg
        # Vectorized SkyCoord
        self.coords = SkyCoord(ra=ras, dec=decs, frame='icrs')

    def get_events_by_altitude(self, altitudes_deg: List[float], start_time: Optional[Union[Time, datetime]] = None) -> \
    Dict[
        str, List[Dict]]:
        t = _ensure_time(start_time)

        times_grid = self._get_time_grid(t)
        frame = AltAz(obstime=times_grid, location=self.location)

        full_results = {}

        # Iterate over stars.
        for idx, star_id in enumerate(self.ids):
            single_star = self.coords[idx]
            star_altaz_over_time = single_star.transform_to(frame)

            star_events = []
            for alt_target in altitudes_deg:
                event_times = self._find_altitude_crossings(alt_target, times_grid, star_altaz_over_time.alt)

                if len(event_times) > 0:
                    frame_exact = AltAz(obstime=event_times, location=self.location)
                    exact_coo = single_star.transform_to(frame_exact)

                    for i, t in enumerate(event_times):
                        star_events.append({
                            'target_alt': alt_target,
                            # FORCE UTC DATETIME
                            'time_utc': t.to_datetime(timezone.utc),
                            'actual_alt': exact_coo.alt.deg[i],
                            'az': exact_coo.az.deg[i],
                            'ra': single_star.ra.deg,
                            'dec': single_star.dec.deg,
                            'body': star_id
                        })
            star_events.sort(key=lambda x: x['time_utc'])
            full_results[star_id] = star_events

        return full_results

    def get_ephemeris(self, times: Union[List[Time], List[datetime], Time, datetime]) -> Dict[str, List[Dict]]:
        t = _ensure_time(times, wrap_scalar=True)

        results = {sid: [] for sid in self.ids}

        # Iterate over time steps
        for ti in t:
            frame_i = AltAz(obstime=ti, location=self.location)
            # Transform all stars for this specific time
            current_altaz = self.coords.transform_to(frame_i)

            for idx, sid in enumerate(self.ids):
                results[sid].append({
                    # FORCE UTC DATETIME
                    'time_utc': ti.to_datetime(timezone.utc),
                    'alt': current_altaz[idx].alt.deg,
                    'az': current_altaz[idx].az.deg,
                    'ra': self.coords[idx].ra.deg,
                    'dec': self.coords[idx].dec.deg,
                    'body': sid
                })

        return results


# ==========================================
# Convenience Functions
# ==========================================


def calculate_sun_rise_set(date: datetime, horiz_height: float, sunrise: bool,
                           latitude: float, longitude: float,
                           elevation: float) -> Optional[datetime]:
    """Next sunrise or sunset at the given horizon altitude.

    Convenience wrapper around :py:meth:`Sun.get_next_event_by_altitude`.
    Returning the next crossing in a chosen direction is a special case
    of "any altitude event of the Sun", so this delegates to the
    general primitive instead of calling astroplan directly. That keeps
    the algorithm consistent with the rest of pyaraucaria's ephemeris
    (5-minute grid + cubic-spline root finding) and — crucially —
    guarantees a real ``datetime`` (or ``None``) return type rather
    than the ``MaskedNDArray`` astroplan emits when no crossing exists
    in the next 24 h.

    :param date: UTC start time of the search.
    :param horiz_height: altitude in degrees (0 for the visible horizon,
        -18 for astronomical twilight, etc.).
    :param sunrise: ``True`` for the next time the sun rises through
        ``horiz_height``; ``False`` for the next time it sets through
        it.
    :param latitude: observer latitude in degrees.
    :param longitude: observer longitude in degrees.
    :param elevation: observer elevation in metres.
    :return: UTC ``datetime`` of the next event, or ``None`` if no
        such event occurs within 24 h of ``date`` (e.g. polar geometry,
        or an altitude the sun never reaches at this site/season).
    """
    location = EarthLocation(lat=latitude * u.deg,
                             lon=longitude * u.deg,
                             height=elevation * u.m)
    direction = 'rising' if sunrise else 'setting'
    return Sun(location).get_next_event_by_altitude(
        horiz_height, direction, start_time=date)


def calculate_moon_rise_set(date: datetime, horiz_height: float, moonrise: bool,
                            latitude: float, longitude: float,
                            elevation: float) -> Optional[datetime]:
    """Next moonrise or moonset at the given horizon altitude.

    Mirror of :py:func:`calculate_sun_rise_set` for the Moon — same
    contract: returns a UTC ``datetime`` for the next crossing in the
    chosen direction, or ``None`` if no such crossing occurs within 24 h
    (e.g. the moon stays below the requested altitude all night).
    Delegates to :py:meth:`Moon.get_next_event_by_altitude` so both
    convenience wrappers share one algorithm.

    :param date: UTC start time of the search.
    :param horiz_height: altitude in degrees (0 for the visible horizon).
    :param moonrise: ``True`` for the next time the moon rises through
        ``horiz_height``; ``False`` for the next time it sets through it.
    :param latitude: observer latitude in degrees.
    :param longitude: observer longitude in degrees.
    :param elevation: observer elevation in metres.
    :return: UTC ``datetime`` of the next event, or ``None`` if no such
        event occurs within 24 h of ``date``.
    """
    location = EarthLocation(lat=latitude * u.deg,
                             lon=longitude * u.deg,
                             height=elevation * u.m)
    direction = 'rising' if moonrise else 'setting'
    return Moon(location).get_next_event_by_altitude(
        horiz_height, direction, start_time=date)


def moon_separation(ra: float, dec: float, utc_time: Time):
    """
    Func. returns moon separation
    :param ra: RA in float (J2000 / ICRS)
    :param dec: Dec in float (J2000 / ICRS)
    :param utc_time: time of calculation in astropy Time format
    :return: Separation in deg
    """
    moon = get_body(body="moon", time=utc_time)
    obj_coo = SkyCoord(ra=ra * u.deg, dec=dec * u.deg, frame='icrs') # FK5(equinox=utc_time))
    return float(moon.separation(obj_coo).to(u.deg).deg)


def moon_phase(date_utc: datetime, latitude: float = None, longitude: float = None, elevation: float = None):
    """
    Return moon phase (illumination) in percent (0–100).

    Moon illumination is a geocentric quantity; *latitude*, *longitude*, and
    *elevation* are therefore not used in the calculation.  Passing them is
    deprecated and they will be removed in a future release.

    :param date_utc: UTC date/time of the calculation
    :param latitude: **deprecated** – ignored, kept for backward compatibility
    :param longitude: **deprecated** – ignored, kept for backward compatibility
    :param elevation: **deprecated** – ignored, kept for backward compatibility
    :return: moon phase (illumination) in % (range 0–100)
    """
    if latitude is not None or longitude is not None or elevation is not None:
        warnings.warn(
            "latitude, longitude, and elevation are deprecated in moon_phase() "
            "and are ignored; moon illumination is a geocentric quantity. "
            "Pass only date_utc in new code.",
            DeprecationWarning,
            stacklevel=2,
        )
    return float(moon_illumination(_ensure_time(date_utc))) * 100
