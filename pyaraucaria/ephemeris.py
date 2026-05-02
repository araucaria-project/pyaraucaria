from typing import List, Dict, Union, Optional
import sys
import csv
import warnings
from datetime import datetime, timezone
import argparse

import numpy as np
from scipy.interpolate import UnivariateSpline

from astropy import units as u
from astropy.time import Time
from astropy.coordinates import SkyCoord, EarthLocation, AltAz, get_body, FK5
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

    def get_next_event_by_altitude(self, altitude_deg: float, direction: str,
                                   start_time: Optional[Union[Time, datetime]] = None
                                   ) -> Optional[datetime]:
        """First crossing of ``altitude_deg`` after ``start_time`` going
        in the requested vertical direction.

        Generic primitive built on :py:meth:`get_events_by_altitude` —
        the convenience function :py:func:`calculate_sun_rise_set` is a
        thin wrapper around this for the Sun, and any caller that wants
        e.g. "next moonset at horizon" should use this rather than
        re-implementing the rise/set selection.

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
        events = self.get_events_by_altitude([altitude_deg], start_time=start_time)
        if not events:
            return None
        # `get_events_by_altitude` returns crossings without direction
        # info. We derive it from the body's altitude at start_time:
        # each crossing flips the above/below-altitude state, so the
        # first crossing after starting "above" is a setting, after
        # starting "below" is a rising — and so on alternately.
        start = _ensure_time(start_time)
        alt_now = float(self.get_ephemeris(start)[0]['alt'])
        state_above = alt_now > altitude_deg
        want_rising = (direction == 'rising')
        for ev in events:
            is_rising = not state_above
            if is_rising == want_rising:
                return ev['time_utc']
            state_above = not state_above
        return None


class Sun(CelestialBody):
    """
    Class for Solar calculations.
    """

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
