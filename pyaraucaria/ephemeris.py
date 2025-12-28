from typing import List, Dict, Union, Optional
import sys
import csv
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


class Sun(CelestialBody):
    """
    Class for Solar calculations.
    """

    def get_events_by_altitude(self, altitudes_deg: List[float], start_time: Optional[Union[Time, datetime]] = None) -> \
    List[Dict]:
        if start_time is None:
            start_time = Time.now()
        else:
            start_time = Time(start_time)

        # Generate coarse grid for the next 24 hours
        times_grid = self._get_time_grid(start_time)
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
        t = Time(times)
        # Handle single scalar time vs list of times
        is_scalar = t.isscalar
        if is_scalar:
            t = Time([t])

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
        if start_time is None:
            start_time = Time.now()
        else:
            start_time = Time(start_time)

        times_grid = self._get_time_grid(start_time)
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
        t = Time(times)
        is_scalar = t.isscalar
        if is_scalar:
            t = Time([t])

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
        t = Time(times)
        is_scalar = t.isscalar
        if is_scalar:
            t = Time([t])

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
        if start_time is None:
            start_time = Time.now()
        else:
            start_time = Time(start_time)

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
        t = Time(times)
        is_scalar = t.isscalar
        if is_scalar:
            t = Time([t])

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
        if start_time is None:
            start_time = Time.now()
        else:
            start_time = Time(start_time)

        times_grid = self._get_time_grid(start_time)
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
        t = Time(times)
        is_scalar = t.isscalar
        if is_scalar:
            t = Time([t])

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
                           latitude: float, longitude: float, elevation: float):
    """
    Calculate next sunrise or sunset at horizon height
    :param date: utc date of start calculating
    :param horiz_height: the height over (or under) horizon
    :param sunrise: if true the sunrise will be calculated, else sunset
    :param latitude: latitude of observer
    :param longitude: longitude of observer
    :param elevation: elevation of observer
    :return: utc datetime time of next sunrise / sunset
    """
    date = Time(val=date)
    obs = Observer(latitude=latitude,
                   longitude=longitude,
                   elevation=elevation * u.m)
    if sunrise:
        return obs.sun_rise_time(date, which='next', horizon=horiz_height * u.deg).to_datetime(timezone.utc)
    else:
        return obs.sun_set_time(date, which='next', horizon=horiz_height * u.deg).to_datetime(timezone.utc)


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


def moon_phase(date_utc: datetime, latitude: float, longitude: float, elevation: float):
    """
    Func. returns moon phase (illumination).
    :param date_utc: calculating utc date
    :param latitude: latitude of observer
    :param longitude: longitude of observer
    :param elevation: elevation of observer
    :return: moon phase (illumination) in % (range 0-100)
    """

    date = Time(val=date_utc)
    obs = Observer(latitude=latitude,
                   longitude=longitude,
                   elevation=elevation * u.m)
    return obs.moon_illumination(time=date) * 100
