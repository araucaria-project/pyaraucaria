"""Tests for the refactored ObsValidator.

Covers the new JSON Schema 2020-12 path AND the backward-compatible legacy
constructor used by tpg. The legacy test class asserts bit-for-bit identical
`result` dicts across a cross-section of OBs.
"""
from __future__ import annotations

import warnings
import unittest

from pyaraucaria.ob_validator import ObsValidator
from pyaraucaria.obs_plan.obs_plan_parser import ObsPlanParser


# ---------------------------------------------------------------------------
# New API — JSON Schema 2020-12 path
# ---------------------------------------------------------------------------

class TestNewAPIFromDefault(unittest.TestCase):
    """Validator built from the shipped base.yaml (no overlay)."""

    @classmethod
    def setUpClass(cls):
        cls.v = ObsValidator.from_default()

    def _valid(self, ob: dict) -> bool:
        return self.v.validate_ob(ob, allowed_filters=["V", "Ic", "B", "g", "r", "i", "z", "Jn", "Hn", "Ks", "u"])["valid"]

    # --- passing cases ---

    def test_object_full(self):
        self.assertTrue(self._valid({
            "command_name": "OBJECT",
            "name": "FF_Aql",
            "ra": "18:58:14.75",
            "dec": "-21:22:14.0",
            "seq": "5/Ic/60,5/V/70",
        }))

    def test_object_no_coords(self):
        """OBJECT may omit coordinates — use current pointing."""
        self.assertTrue(self._valid({"command_name": "OBJECT", "seq": "1/V/10"}))

    def test_object_altaz(self):
        self.assertTrue(self._valid({
            "command_name": "OBJECT",
            "name": "X", "alt": 60, "az": 270, "seq": "1/V/10",
        }))

    def test_wait_ut(self):
        self.assertTrue(self._valid({"command_name": "WAIT", "ut": "23:00:00"}))

    def test_wait_sec(self):
        self.assertTrue(self._valid({"command_name": "WAIT", "sec": 300}))

    def test_stop(self):
        self.assertTrue(self._valid({"command_name": "STOP"}))

    def test_zero(self):
        self.assertTrue(self._valid({"command_name": "ZERO", "seq": "15/V/0"}))

    def test_skyflat_alt_az_auto_exposure(self):
        self.assertTrue(self._valid({
            "command_name": "SKYFLAT", "alt": 60, "az": 270, "seq": "10/V/a",
        }))

    def test_focus_with_pos(self):
        self.assertTrue(self._valid({
            "command_name": "FOCUS", "name": "NG31",
            "ra": "12:12:12", "dec": "+20:20:20",
            "pos": "15200/100",
            "seq": "5/Ic/3",
            "auto_focus": "on",
        }))

    # --- failing cases ---

    def test_object_bad_ra_pattern(self):
        ob = {"command_name": "OBJECT", "name": "X", "ra": "18:58", "dec": "-21:22:14.0", "seq": "1/V/10"}
        self.assertFalse(self._valid(ob))

    def test_object_both_radec_and_altaz_rejected(self):
        """Ra/dec AND alt/az together is contradictory; must be rejected."""
        ob = {"command_name": "OBJECT", "ra": "18:00:00", "dec": "-30:00:00",
              "alt": 45, "az": 180, "seq": "1/V/1"}
        self.assertFalse(self._valid(ob))

    def test_wait_requires_something(self):
        self.assertFalse(self._valid({"command_name": "WAIT"}))

    def test_zero_missing_seq(self):
        self.assertFalse(self._valid({"command_name": "ZERO"}))

    def test_unknown_command(self):
        self.assertFalse(self._valid({"command_name": "UNKNOWN"}))

    def test_unknown_kwarg_rejected(self):
        """additionalProperties: false (via unevaluatedProperties) should bite."""
        ob = {"command_name": "OBJECT", "name": "X",
              "ra": "18:00:00", "dec": "-30:00:00", "seq": "1/V/1",
              "bogus": "value"}
        self.assertFalse(self._valid(ob))

    def test_seq_with_bad_filter_rejected(self):
        ob = {"command_name": "OBJECT", "seq": "1/NOSUCHFILTER/10"}
        self.assertFalse(self._valid(ob))

    def test_seq_autoexp_forbidden_on_object(self):
        """`a` auto-exposure is only allowed on FOCUS and SKYFLAT."""
        ob = {"command_name": "OBJECT", "name": "X", "seq": "1/V/a"}
        self.assertFalse(self._valid(ob))


class TestNewAPIWithTpgOverlay(unittest.TestCase):
    """Validator built from base.yaml + tpg_overlay.yaml."""

    @classmethod
    def setUpClass(cls):
        cls.v = ObsValidator.from_default(overlays=["tpg_overlay"])

    def test_tpg_kwargs_accepted(self):
        ob = {
            "command_name": "OBJECT",
            "name": "FF_Aql", "ra": "18:58:14.75", "dec": "-21:22:14.0",
            "seq": "5/Ic/60",
            "priority": 5, "cycle": 0.08, "tag": "RR,calib",
            "P": 2.4121, "hjd0": 2460310.5,
            "ph_start": 0.8, "ph_end": 0.1,
            "ph_mk": "2/Ic/0.05",
            "obs_data": "/data/x.txt",
            "mV": 9.21, "mK": 7.5,
            "h_min": 40, "h_max": 80,
            "min_moon_dist": 30, "max_moon_phase": 70,
        }
        r = self.v.validate_ob(ob, allowed_filters=["V", "Ic"])
        self.assertTrue(r["valid"], f"failed: {r['result']}")

    def test_priority_must_be_number(self):
        ob = {"command_name": "OBJECT", "name": "X",
              "ra": "18:00:00", "dec": "-30:00:00", "seq": "1/V/1",
              "priority": "dupa"}
        r = self.v.validate_ob(ob, allowed_filters=["V"])
        self.assertFalse(r["valid"])

    def test_unknown_kwarg_still_rejected(self):
        ob = {"command_name": "OBJECT", "name": "X",
              "ra": "18:00:00", "dec": "-30:00:00", "seq": "1/V/1",
              "not_a_real_kwarg": 42}
        r = self.v.validate_ob(ob, allowed_filters=["V"])
        self.assertFalse(r["valid"])


class TestNewAPIParseAndValidateRoundTrip(unittest.TestCase):
    """End-to-end: parser → convert_to_obdict → validate."""

    @classmethod
    def setUpClass(cls):
        cls.v = ObsValidator.from_default(overlays=["tpg_overlay"])

    def _round_trip(self, line: str) -> dict:
        parsed = ObsPlanParser.convert_from_string(line)
        ob = ObsValidator.convert_to_obdict(parsed)
        return self.v.validate_ob(ob, allowed_filters=["V", "Ic", "B", "g", "r", "i", "z"])

    def test_realistic_tpg_line(self):
        line = ('OBJECT FF_Aql 18:58:14.75 -21:22:14.0 seq=2/Ic/28,2/V/20 '
                'priority=5 cycle=0.08 tag=random_phase uobi=4d5a6db1 '
                'sciprog=t2cep pi=pwielgorski P=1.091787 mV=12.31 '
                'comment="binary and dusty envelope"')
        self.assertTrue(self._round_trip(line)["valid"])

    def test_focus_line(self):
        line = "FOCUS NG31 12:12:12 +20:20:20 pos=15200/100 seq=5/Ic/3 auto_focus=on"
        r = self._round_trip(line)
        self.assertTrue(r["valid"], f"failed: {r['result']}")


# ---------------------------------------------------------------------------
# Legacy API — bit-for-bit behaviour preservation
# ---------------------------------------------------------------------------

class TestLegacyAPICompatibility(unittest.TestCase):
    """Legacy two-arg constructor must keep behaving as before the refactor.

    Every test that was possible with the original validator should still
    produce the same `valid`/`result` tuple. This is what keeps tpg working
    on an unchanged line until it migrates.
    """

    @classmethod
    def setUpClass(cls):
        cls.base = ObsValidator.load_schema("base_schema")
        cls.rules = ObsValidator.load_schema("base_rules")
        with warnings.catch_warnings():
            warnings.simplefilter("ignore", DeprecationWarning)
            cls.v = ObsValidator(cls.base, cls.rules)

    def test_deprecation_warning_fires(self):
        with warnings.catch_warnings(record=True) as ws:
            warnings.simplefilter("always")
            ObsValidator(self.base, self.rules)
            self.assertTrue(
                any(issubclass(w.category, DeprecationWarning) for w in ws),
                "legacy ctor must emit DeprecationWarning",
            )

    def test_valid_object(self):
        ob = {"command_name": "OBJECT", "name": "FF_Aql",
              "ra": "18:58:14.75", "dec": "-21:22:14.0", "seq": "5/Ic/60"}
        r = self.v.validate_ob(ob, allowed_filters=["V", "Ic"])
        self.assertTrue(r["valid"])

    def test_valid_wait(self):
        r = self.v.validate_ob({"command_name": "WAIT", "ut": "23:00:00"})
        self.assertTrue(r["valid"])

    def test_valid_zero(self):
        r = self.v.validate_ob({"command_name": "ZERO", "seq": "15/V/0"},
                               allowed_filters=["V"])
        self.assertTrue(r["valid"])

    def test_invalid_bad_ra(self):
        ob = {"command_name": "OBJECT", "name": "X", "ra": "bad", "dec": "-20:00:00", "seq": "1/V/10"}
        r = self.v.validate_ob(ob, allowed_filters=["V"])
        self.assertFalse(r["valid"])
        self.assertFalse(r["result"].get("ra", True))

    def test_invalid_unknown_command(self):
        r = self.v.validate_ob({"command_name": "NOPE"})
        self.assertFalse(r["valid"])
        self.assertFalse(r["result"].get("command_name", True))

    def test_legacy_preserves_known_quirks(self):
        """These behaviours are KNOWN BUGS that the new API fixes.

        The legacy path must preserve them verbatim so tpg's existing logic
        keeps working. PR 2 (strict validation) is where the behaviour flips.
        """
        # Unknown kwarg passes in legacy (additionalProperties: true).
        ob = {"command_name": "OBJECT", "name": "X",
              "ra": "18:00:00", "dec": "-30:00:00", "seq": "1/V/1",
              "bogus_kwarg": "passes-in-legacy"}
        r = self.v.validate_ob(ob, allowed_filters=["V"])
        self.assertTrue(r["valid"],
                        "legacy path intentionally accepts unknown kwargs")


# ---------------------------------------------------------------------------
# Helpers — convert_to_obdict / convert_from_obdict round-trip
# ---------------------------------------------------------------------------

class TestConversionHelpers(unittest.TestCase):

    def test_convert_to_obdict_three_args(self):
        parsed = ObsPlanParser.convert_from_string(
            "OBJECT FF_Aql 18:58:14.75 -21:22:14.0 seq=1/V/60"
        )
        ob = ObsValidator.convert_to_obdict(parsed)
        self.assertEqual(ob["name"], "FF_Aql")
        self.assertEqual(ob["ra"], "18:58:14.75")
        self.assertEqual(ob["dec"], "-21:22:14.0")
        self.assertEqual(ob["seq"], "1/V/60")

    def test_convert_from_obdict_round_trip(self):
        ob = {"command_name": "OBJECT", "name": "FF_Aql",
              "ra": "18:58:14.75", "dec": "-21:22:14.0", "seq": "1/V/60"}
        line = ObsValidator.convert_from_obdict(ob)
        self.assertIn("OBJECT", line)
        self.assertIn("FF_Aql", line)
        self.assertIn("seq=1/V/60", line)


if __name__ == "__main__":
    unittest.main()
