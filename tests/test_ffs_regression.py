"""Regression (snapshot) tests for ``pyaraucaria.ffs``.

These tests do NOT check that FFS is right -- they check that it did not
change. Each case replays a call sequence used in production (TOI, Focus,
FitsView) on deterministic synthetic frames and compares every output with
the values stored in ``tests/data/ffs_snapshot.json``, known bugs included.

A failure here means production results would change. If the change is
intended, regenerate the snapshot and commit it together with the code:

    python -m tests.test_ffs_regression --regen
"""

import hashlib
import json
import os
import sys
import unittest
import warnings

import numpy as np

from pyaraucaria.ffs import FFS
from tests.ffs_synth import SCENARIOS, scenario

SNAPSHOT = os.path.join(os.path.dirname(__file__), "data", "ffs_snapshot.json")
RTOL = 1e-7
ATOL = 1e-9

STAR_COLUMNS = ["x", "y", "max_adu", "box_mag", "bkg", "fwhm", "fwhm_x", "fwhm_y",
                "ellipticity", "theta", "cpe", "shape", "ci"]


def _mask_digest(mask):
    mask = np.asarray(mask, dtype=bool)
    return {"count": int(mask.sum()),
            "sha256": hashlib.sha256(np.packbits(mask).tobytes()).hexdigest()}


def _stars(ffs):
    return {c: ffs.stars[c] for c in STAR_COLUMNS if c in ffs.stars.colnames}


# --- Production call sequences ---------------------------------------------

def run_toi(image):
    # TOI/fits_gui.py FFS_Worker.run
    ffs = FFS(image)
    ffs.saturation = 45000
    ffs.calc_frame_fwhm(threshold=10, fwhm=5, box=15, N_stars=50, clip=4)
    ffs.sky_gradient(n_segments=7)
    return {"frame": ffs.stats["frame"], "stars": _stars(ffs), "sky": ffs.stats["sky"]}


def run_focus(image):
    # pyaraucaria/focus.py Focus.fwhm, called with its defaults
    from pyaraucaria.focus import Focus
    frame_fwhm = Focus.fwhm(image, saturation=45000)
    return {"frame_fwhm": np.nan if frame_fwhm is None else frame_fwhm}


def run_fitsview(image):
    # FitsView/fitsview/FitsView_widgets.py FFSWindow (init, find, sky, lines)
    ffs = FFS(image)
    ffs.saturation = 45000
    ffs.mk_stats()
    ffs.find_stars(threshold=5, fwhm=4)
    ffs.calc_frame_fwhm(threshold=5, fwhm=4, box=8, N_stars=20)
    ffs.star_info(box=8, N_stars=None)
    out = {"frame": dict(ffs.stats["frame"]), "stars": _stars(ffs)}
    ffs.sky_gradient(n_segments=10)
    out["sky"] = ffs.stats["sky"]
    ffs.find_lines()
    out["lines"] = ffs.stats["lines"]
    out["maska"] = _mask_digest(ffs.masks["threshold"])
    return out


def run_find_stars_default(image):
    ffs = FFS(image)
    ffs.mk_stats()
    ffs.find_stars()
    return {"coo": ffs.coo, "adu": ffs.adu,
            "sigma": ffs.fs_sigma, "kernel_sigma": ffs.fs_kernel_sigma}


def run_find_stars_filters(image):
    # The detection options added in 82a233e (guider use).
    out = {}
    for rank_by in ("raw", "smoothed", "aperture"):
        ffs = FFS(image)
        ffs.mk_stats()
        ffs.find_stars(threshold=5, fwhm=4, min_smoothed_sigma=5, rank_by=rank_by,
                       max_concentration=0.5, min_pixels_above_threshold=5)
        out[rank_by] = {"coo": ffs.coo, "adu": ffs.adu}
    return out


def run_line_filter(image):
    return {"mask": _mask_digest(FFS.line_filter(image.astype(float)))}


CASES = {
    "toi": run_toi,
    "focus": run_focus,
    "fitsview": run_fitsview,
    "find_stars_default": run_find_stars_default,
    "find_stars_filters": run_find_stars_filters,
    "line_filter": run_line_filter,
}


def compute(case, scenario_name):
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", RuntimeWarning)
        return _plain(CASES[case](scenario(scenario_name).image))


def _plain(obj):
    """Convert results to JSON-compatible python types."""
    if isinstance(obj, dict):
        return {str(k): _plain(v) for k, v in obj.items()}
    if isinstance(obj, (list, tuple)):
        return [_plain(v) for v in obj]
    if hasattr(obj, "tolist"):          # numpy arrays, scalars, astropy Column
        return _plain(obj.tolist()) if np.ndim(obj) else obj.tolist()
    return obj


# --- Comparison -------------------------------------------------------------

def assert_same(test, expected, actual, path=""):
    if isinstance(expected, dict):
        test.assertIsInstance(actual, dict, path)
        test.assertEqual(sorted(expected), sorted(actual), f"keys differ at {path or '/'}")
        for k in expected:
            assert_same(test, expected[k], actual[k], f"{path}/{k}")
    elif isinstance(expected, str):
        test.assertEqual(expected, actual, path)
    elif isinstance(expected, list) and any(isinstance(v, (dict, str)) for v in expected):
        test.assertEqual(len(expected), len(actual), f"length differs at {path}")
        for i, (e, a) in enumerate(zip(expected, actual)):
            assert_same(test, e, a, f"{path}[{i}]")
    else:
        e = np.asarray(expected, dtype=float)
        a = np.asarray(actual, dtype=float)
        test.assertEqual(e.shape, a.shape, f"shape differs at {path}")
        np.testing.assert_allclose(a, e, rtol=RTOL, atol=ATOL, equal_nan=True,
                                   err_msg=f"values differ at {path}")


class TestFFSRegression(unittest.TestCase):

    @classmethod
    def setUpClass(cls):
        if not os.path.exists(SNAPSHOT):
            raise unittest.SkipTest(
                f"{SNAPSHOT} missing; run: python -m tests.test_ffs_regression --regen")
        with open(SNAPSHOT) as fh:
            cls.snapshot = json.load(fh)

    def test_snapshot_covers_all_cases(self):
        expected = {f"{c}:{s}" for c in CASES for s in SCENARIOS}
        self.assertEqual(expected, set(self.snapshot["results"]))

    def _check(self, case):
        for name in SCENARIOS:
            with self.subTest(scenario=name):
                expected = self.snapshot["results"][f"{case}:{name}"]
                assert_same(self, expected, compute(case, name), f"{case}:{name}")

    def test_toi(self):
        self._check("toi")

    def test_focus(self):
        self._check("focus")

    def test_fitsview(self):
        self._check("fitsview")

    def test_find_stars_default(self):
        self._check("find_stars_default")

    def test_find_stars_filters(self):
        self._check("find_stars_filters")

    def test_line_filter(self):
        self._check("line_filter")


def regenerate():
    results = {f"{c}:{s}": compute(c, s) for c in CASES for s in SCENARIOS}
    data = {
        "note": "Generated by `python -m tests.test_ffs_regression --regen`. "
                "Records current FFS behaviour (bugs included); see module docstring.",
        "numpy": np.__version__,
        "results": results,
    }
    with open(SNAPSHOT, "w") as fh:
        json.dump(data, fh, indent=1, sort_keys=True)
    print(f"wrote {SNAPSHOT} ({len(results)} cases)")


if __name__ == "__main__":
    if "--regen" in sys.argv:
        regenerate()
    else:
        unittest.main()
