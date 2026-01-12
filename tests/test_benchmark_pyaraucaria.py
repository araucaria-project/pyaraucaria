import unittest
import time
from datetime import datetime, timezone
from astropy.coordinates import EarthLocation
from astropy.time import Time
import astropy.units as u

# Reference Implementation
# Ensure these imports match your project structure
from pyaraucaria.ephemeris import calculate_sun_rise_set, moon_phase, Sun, Moon


class TestBenchmarks(unittest.TestCase):

    def setUp(self):
        """
        Set up common test fixtures (equivalent to pytest fixtures).
        """
        # OCA approx coordinates
        self.lat = -24.6
        self.lon = -70.2
        self.elev = 2800.0

        # EarthLocation object for ocacal
        self.earth_location = EarthLocation(
            lat=self.lat * u.deg,
            lon=self.lon * u.deg,
            height=self.elev * u.m
        )

        # A fixed datetime for comparison
        self.check_time = datetime(2025, 5, 20, 12, 0, 0, tzinfo=timezone.utc)

    def test_benchmark_sun_events(self):
        """
        Compare calculation of Sunset/Sunrise times and execution speed.
        """
        print("\n\n--- BENCHMARK: Sun Rise/Set ---")
        print(f"check_time: {self.check_time.isoformat()} UTC")

        # 1. Pyaraucaria Execution
        # ------------------------
        t0 = time.perf_counter()

        # Calculate Sunrise (horizon=0)
        p_rise = calculate_sun_rise_set(
            self.check_time, 0.0, True, self.lat, self.lon, self.elev
        )
        # Calculate Sunset (horizon=0)
        p_set = calculate_sun_rise_set(
            self.check_time, 0.0, False, self.lat, self.lon, self.elev
        )

        t_py = (time.perf_counter() - t0) * 1000  # ms
        print(f"Pyaraucaria time: {t_py:.3f} ms")
        print(f"Pyaraucaria Rise: {p_rise}")
        print(f"Pyaraucaria Set : {p_set}")

        # 2. Ocacal Execution
        # -------------------
        t0 = time.perf_counter()

        sun = Sun(self.earth_location)
        # We ask for 0 degrees altitude crossing
        events = sun.get_events_by_altitude([0.0], start_time=self.check_time)

        t_oca = (time.perf_counter() - t0) * 1000  # ms
        print(f"Ocacal time     : {t_oca:.3f} ms")

        # 3. Validation & Comparison
        # --------------------------

        # Filter Ocacal results
        # We assume:
        # - If Azimuth is East (< 180), it's rising.
        # - If Azimuth is West (> 180), it's setting.
        o_rise = next((e for e in events if e['az'] < 180), None)
        o_set = next((e for e in events if e['az'] > 180), None)

        print(f"Ocacal Rise     : {o_rise['time_utc'] if o_rise else 'Not found'}")
        print(f"Ocacal Set      : {o_set['time_utc'] if o_set else 'Not found'}")

        self.assertIsNotNone(o_rise, "Ocacal failed to find sunrise")
        self.assertIsNotNone(o_set, "Ocacal failed to find sunset")

        # Calculate difference in seconds
        diff_rise = abs((o_rise['time_utc'] - p_rise).total_seconds())
        diff_set = abs((o_set['time_utc'] - p_set).total_seconds())

        print(f"Delta Rise      : {diff_rise:.2f} seconds")
        print(f"Delta Set       : {diff_set:.2f} seconds")

        # Assertions
        # We allow a small tolerance (e.g., 60 seconds) because implementations differ:
        # - Ocacal uses cubic spline interpolation on a 5-min grid.
        # - Pyaraucaria uses astroplan's built-in grid search/interpolation.
        self.assertLess(diff_rise, 60.0, f"Sunrise differs by {diff_rise}s")
        self.assertLess(diff_set, 60.0, f"Sunset differs by {diff_set}s")

        # Optional: Speed comparison warning (not failure)
        if t_oca > t_py:
            print(f"-> Note: Ocacal was {t_oca / t_py:.1f}x slower (expected due to spline setup overhead)")
        else:
            print(f"-> Note: Ocacal was faster!")

    def test_benchmark_moon_phase(self):
        """
        Compare calculation of Moon Phase and execution speed.
        """
        print("\n\n--- BENCHMARK: Moon Phase ---")

        # 1. Pyaraucaria Execution
        # ------------------------
        t0 = time.perf_counter()

        p_phase = moon_phase(self.check_time, self.lat, self.lon, self.elev) * 0.01

        t_py = (time.perf_counter() - t0) * 1000  # ms
        print(f"Pyaraucaria time : {t_py:.3f} ms")
        print(f"Pyaraucaria Phase: {p_phase:.6f}")

        # 2. Ocacal Execution
        # -------------------
        t0 = time.perf_counter()

        moon = Moon(self.earth_location)
        # get_ephemeris takes a list or scalar Time
        # result = moon.get_ephemeris(Time(self.check_time))
        result = moon.get_phase(self.check_time)
        o_phase = result[0]['phase']

        t_oca = (time.perf_counter() - t0) * 1000  # ms
        print(f"Ocacal time      : {t_oca:.3f} ms")
        print(f"Ocacal Phase     : {o_phase:.6f}")

        # 3. Validation & Comparison
        # --------------------------
        diff = abs(p_phase - o_phase)
        print(f"Delta Phase      : {diff:.6f}")

        # Assertions
        # Should be nearly identical as both likely wrap astroplan.moon.moon_illumination
        # or similar astropy functions.
        self.assertLess(diff, 0.001, f"Moon phase differs by {diff}")

        if t_oca > t_py:
            print(f"-> Note: Ocacal was {t_oca / t_py:.1f}x slower")
        else:
            print(f"-> Note: Ocacal was faster!")


if __name__ == '__main__':
    unittest.main()