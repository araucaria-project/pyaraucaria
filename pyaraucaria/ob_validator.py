"""Observing-block validator.

Two usage styles are supported:

**New API (preferred).**  Pass a single JSON Schema 2020-12 document; `$ref`
across files is resolved from `pyaraucaria/schemas/`. Validation uses
`jsonschema.Draft202012Validator` directly.

    from pyaraucaria.ob_validator import ObsValidator
    v = ObsValidator.from_default()                         # base.yaml
    v = ObsValidator.from_default(overlays=["tpg_overlay"]) # base + tpg extras
    v = ObsValidator(my_schema_dict)                        # custom schema
    r = v.validate_ob({"command_name": "OBJECT", ...})

**Legacy API (deprecated — kept for backward compatibility with tpg).**  The
two-arg constructor takes `(base_schema, command_rules)` as loaded from the
old `base_schema.yaml` + `base_rules.yaml` pair. Behaviour is preserved
bit-for-bit from the pre-refactor implementation. Emits a `DeprecationWarning`
on construction.

    v = ObsValidator(BASE_SCHEMA, COMMAND_RULES)
    r = v.validate_ob(ob)
"""
from __future__ import annotations

import os
import re
import warnings
from copy import deepcopy
from functools import lru_cache
from pathlib import Path
from typing import Any, Iterable

import yaml
from jsonschema import Draft7Validator, Draft202012Validator
from referencing import Registry, Resource
from referencing.jsonschema import DRAFT202012

from pyaraucaria.obs_plan.seq import validate_seq as _validate_seq


_SCHEMAS_DIR = Path(__file__).parent / "schemas"


# ---------------------------------------------------------------------------
# Schema loading / $ref resolution
# ---------------------------------------------------------------------------

@lru_cache(maxsize=None)
def _load_yaml_cached(abs_path: str) -> dict:
    with open(abs_path) as f:
        return yaml.safe_load(f)


def _load_yaml(path: Path) -> dict:
    return deepcopy(_load_yaml_cached(str(path.resolve())))


def _build_registry() -> Registry:
    """Populate a referencing Registry with every YAML schema under schemas/.

    Each schema is registered under both its `$id` (when present) and its
    absolute `file://` URI, so `$ref`s written as relative paths resolve
    against the filesystem.
    """
    registry = Registry()
    for p in sorted(_SCHEMAS_DIR.rglob("*.yaml")):
        data = _load_yaml(p)
        resource = Resource(contents=data, specification=DRAFT202012)
        uri = "file://" + str(p.resolve())
        registry = registry.with_resource(uri=uri, resource=resource)
        if isinstance(data, dict) and "$id" in data:
            registry = registry.with_resource(uri=data["$id"], resource=resource)
    return registry


_REGISTRY: Registry | None = None


def _registry() -> Registry:
    global _REGISTRY
    if _REGISTRY is None:
        _REGISTRY = _build_registry()
    return _REGISTRY


def _inline_refs(node: Any, base_dir: Path) -> Any:
    """Resolve relative-path $refs by inlining them.

    Used instead of the Registry when a schema lives outside `schemas/` (e.g.
    a caller-provided dict). For new-style validation we prefer the Registry.
    """
    if isinstance(node, dict):
        ref = node.get("$ref")
        if isinstance(ref, str) and not ref.startswith("#"):
            path, _, pointer = ref.partition("#")
            target = (base_dir / path).resolve()
            sub = _load_yaml(target)
            for part in (p for p in pointer.split("/") if p):
                sub = sub[part]
            return _inline_refs(sub, target.parent)
        return {k: _inline_refs(v, base_dir) for k, v in node.items()}
    if isinstance(node, list):
        return [_inline_refs(x, base_dir) for x in node]
    return node


# ---------------------------------------------------------------------------
# Validator
# ---------------------------------------------------------------------------

class ObsValidator:
    """Schema-driven observing-block validator.

    The class supports two construction styles; see module docstring.
    """

    _DEPRECATION_MSG = (
        "ObsValidator(base_schema, command_rules) is deprecated. "
        "Use ObsValidator.from_default() or ObsValidator(schema) with a "
        "JSON Schema 2020-12 document (see pyaraucaria/schemas/base.yaml)."
    )

    def __init__(
        self,
        schema_or_base: dict,
        command_rules: dict | None = None,
    ) -> None:
        if command_rules is not None:
            warnings.warn(self._DEPRECATION_MSG, DeprecationWarning, stacklevel=2)
            self._mode = "legacy"
            self._base_schema = schema_or_base
            self._command_rules = command_rules
            # Exposed for BC — callers may read these.
            self.base_schema = schema_or_base
            self.command_rules = command_rules
        else:
            self._mode = "new"
            self._schema = schema_or_base
            self._validator = Draft202012Validator(
                schema_or_base,
                registry=_registry(),
            )

    # ------------------------- factory helpers -------------------------

    @classmethod
    def from_default(
        cls,
        overlays: Iterable[str] | None = None,
        *,
        strict: bool = True,
    ) -> "ObsValidator":
        """Build a validator from the shipped `base.yaml` plus optional overlays.

        Parameters
        ----------
        overlays:
            Iterable of overlay schema names to `allOf`-compose with the base.
            Each name is resolved as `pyaraucaria/schemas/<name>.yaml` (the
            `.yaml` suffix is optional). E.g. `["tpg_overlay"]`.
        strict:
            If True (default) wraps the composite with
            `unevaluatedProperties: false` at the top — i.e. any kwarg not
            declared by the base OR by an overlay is rejected.
        """
        branches: list[dict] = [{"$ref": "./base.yaml"}]
        for name in overlays or ():
            if not name.endswith(".yaml"):
                name = name + ".yaml"
            branches.append({"$ref": f"./{name}"})

        composite: dict = (
            {"$ref": "./base.yaml"} if len(branches) == 1 else {"allOf": branches}
        )
        if strict:
            composite = {"allOf": [composite], "unevaluatedProperties": False}

        # Register composite under a synthetic URI so the Registry can resolve it.
        # But since we pass the composite directly to the validator and the
        # Registry already has base.yaml + overlays, $refs will resolve.
        composite = _inline_refs(composite, _SCHEMAS_DIR)
        return cls(composite)

    # ------------------------- schema I/O ------------------------------

    @staticmethod
    def load_schema(name: str) -> dict:
        """Load a YAML schema from `pyaraucaria/schemas/`.

        Accepts bare names (`"base_schema"`), names with extension
        (`"tpg_overlay.yaml"`), and paths (`"commands/OBJECT.yaml"`).
        """
        if not name.endswith(".yaml"):
            name = name + ".yaml"
        path = _SCHEMAS_DIR / name
        if not path.exists():
            raise FileNotFoundError(f"Schema not found: {path}")
        return _load_yaml(path)

    # ------------------------- round-trip helpers ----------------------

    @staticmethod
    def convert_from_obdict(ob: dict) -> str | None:
        if not isinstance(ob, dict):
            return None
        cmd = ob.get("command_name")
        if not cmd:
            return None

        parts: list[str] = [cmd]
        name = ob.get("name")
        ra = ob.get("ra")
        dec = ob.get("dec")
        if name and not (ra or dec):
            parts.append(str(name))
        elif ra and dec and not name:
            parts.extend([str(ra), str(dec)])
        elif name and ra and dec:
            parts.extend([str(name), str(ra), str(dec)])

        skip = {"command_name", "name", "ra", "dec"}
        for k, v in ob.items():
            if k not in skip:
                parts.append(f"{k}={v}")
        return " ".join(parts)

    @staticmethod
    def convert_to_obdict(ob: dict) -> dict | None:
        """Flatten parser output (nested sequences with args/kwargs) into a flat dict.

        Positional args are unpacked by arity — see the module docstring for
        limitations (the stricter PR replaces this with schema-driven unpacking).
        """
        result: dict = {}
        subcommands = ob.get("subcommands", [])
        if isinstance(subcommands, dict):
            subcommands = [subcommands]

        for sub in subcommands:
            if "command_name" in sub:
                result["command_name"] = sub["command_name"]
            if "kwargs" in sub and isinstance(sub["kwargs"], dict):
                result.update(sub["kwargs"])
            if "args" in sub and isinstance(sub["args"], list):
                args = sub["args"]
                if len(args) == 1:
                    result["name"] = args[0]
                elif len(args) == 2:
                    result["ra"] = args[0]
                    result["dec"] = args[1]
                elif len(args) == 3:
                    result["name"] = args[0]
                    result["ra"] = args[1]
                    result["dec"] = args[2]
        return result

    @staticmethod
    def clean_none(obs: dict) -> dict:
        return {k: v for k, v in obs.items() if v is not None}

    # ------------------------- type coercion ---------------------------

    @staticmethod
    def convert_types(obs: dict, schema: dict) -> dict:
        """Coerce string parser outputs to declared schema types.

        Accepts either new-style (with `$ref`s — inspects the top-level schema's
        properties) or legacy-style schemas. Unknown keys pass through untouched.
        """
        properties = _collect_properties(schema)
        converted: dict = {}

        for key, value in obs.items():
            prop = properties.get(key)
            if not prop:
                converted[key] = value
                continue
            expected = prop.get("type")
            if expected is None:
                # could be oneOf(integer, string) etc — try integer then fall through
                converted[key] = _coerce_union(value, prop)
                continue
            if not isinstance(expected, list):
                expected = [expected]

            new_value = value
            ok = False
            for typ in expected:
                try:
                    if typ == "integer":
                        new_value = int(value); ok = True; break
                    if typ == "number":
                        new_value = float(value); ok = True; break
                    if typ == "boolean":
                        if isinstance(value, bool):
                            new_value = value; ok = True; break
                        if isinstance(value, str):
                            low = value.lower()
                            if low in ("true", "1"):
                                new_value = True; ok = True; break
                            if low in ("false", "0"):
                                new_value = False; ok = True; break
                    elif typ == "string":
                        new_value = str(value); ok = True; break
                    elif typ == "null" and value is None:
                        new_value = None; ok = True; break
                except (ValueError, TypeError):
                    pass
            converted[key] = new_value if ok else value
        return converted

    # ------------------------- main entry point ------------------------

    def validate_ob(
        self,
        obs: dict,
        overrides: dict | None = None,
        allowed_filters: Iterable[str] | None = None,
    ) -> dict:
        """Validate an observing block.

        Returns a dict with keys `valid`, `result`, `data`, `required`, `allowed`.
        """
        if self._mode == "legacy":
            return self._validate_ob_legacy(obs, overrides, allowed_filters)
        return self._validate_ob_new(obs, overrides, allowed_filters)

    # ------------------------- new-mode internals ----------------------

    def _validate_ob_new(
        self,
        obs: dict,
        overrides: dict | None,
        allowed_filters: Iterable[str] | None,
    ) -> dict:
        obs_clean = self.clean_none(obs)
        obs_clean = self.convert_types(obs_clean, self._schema)

        result: dict[str, bool] = {key: True for key in obs_clean}

        for error in self._validator.iter_errors(obs_clean):
            if error.path:
                key = error.path[0]
                if isinstance(key, str) and key in result:
                    result[key] = False
            else:
                # Top-level failure (e.g. oneOf branch mismatch, required).
                # Mark command_name so `valid` goes False even if no per-field
                # error was tagged.
                result.setdefault("command_name", True)
                # Use a generic sentinel key so the caller sees the failure.
                result["_schema"] = False

        # Semantic seq check (filter whitelist + auto-exposure policy) — not
        # expressible in JSON Schema alone.
        if "seq" in obs_clean and obs_clean.get("seq"):
            if not _validate_seq(
                obs_clean["seq"],
                command_name=obs_clean.get("command_name"),
                allowed_filters=allowed_filters,
            ):
                result["seq"] = False

        valid = all(v is True for v in result.values())
        properties = _collect_properties(self._schema)
        allowed_fields = sorted(properties.keys())
        return {
            "valid": valid,
            "result": result,
            "data": obs_clean,
            "required": _collect_required(self._schema, obs_clean.get("command_name")),
            "allowed": allowed_fields,
        }

    # ------------------------- legacy-mode internals -------------------

    def _validate_ob_legacy(
        self,
        obs: dict,
        overrides: dict | None,
        allowed_filters: Iterable[str] | None,
    ) -> dict:
        obs_clean = self.clean_none(obs)
        schema = deepcopy(self._base_schema)

        if overrides:
            for k, v in overrides.items():
                properties = schema.setdefault("properties", {})
                if k in properties:
                    properties[k].update(v)
                else:
                    properties[k] = v
        properties = schema.get("properties", {})

        obs_clean = self.convert_types(obs_clean, schema)

        validator = Draft7Validator(schema)
        result = {key: True for key in obs_clean}
        for error in validator.iter_errors(obs_clean):
            if error.validator in ("required", "oneOf"):
                continue
            if error.path:
                key = error.path[0]
                if key in result:
                    result[key] = False

        cmd_type = obs_clean.get("command_name")
        if cmd_type in self._command_rules:
            result.update(self._legacy_validate_rules(obs_clean, self._command_rules[cmd_type]))
        else:
            result["command_name"] = False

        result.update(self._legacy_validate_seq(obs_clean, allowed_filters=allowed_filters))

        required_fields: list[str] = []
        allowed_fields = list(properties.keys())
        if cmd_type in self._command_rules:
            rules = self._command_rules[cmd_type]
            required_fields = rules.get("required", [])
            allowed_fields = list(set(list(properties.keys()) + rules.get("allowed", [])))

        valid = all(v is True for v in result.values())
        return {
            "valid": valid,
            "result": result,
            "data": obs_clean,
            "required": required_fields,
            "allowed": allowed_fields,
        }

    @staticmethod
    def _legacy_validate_rules(obs: dict, rules: dict) -> dict:
        result: dict = {}
        for key in rules.get("required", []):
            if key not in obs:
                result[key] = None
        for subschema in rules.get("one_of", []):
            if not any(k not in obs for k in subschema):
                continue
            satisfied = any(all(k in obs for k in s) for s in rules.get("one_of", []))
            if not satisfied:
                for s in rules.get("one_of", []):
                    for k in s:
                        result[k] = False
        for group in rules.get("one_of_group", []):
            if not any(k in obs for k in group):
                for k in group:
                    result[k] = False
        return result

    @staticmethod
    def _legacy_validate_seq(obs: dict, allowed_filters=None) -> dict:
        # Preserved bit-for-bit from pre-refactor code path.
        result: dict = {}
        allowed_filters = allowed_filters or []
        seq_str = obs.get("seq")
        cmd = obs.get("command_name")
        if seq_str is None:
            return result

        def _validate(seq_str: str) -> bool:
            seq_str = seq_str.strip()
            match = re.match(r"^(\d+)x\((.+)\)$", seq_str)
            if match:
                repeat_count = int(match.group(1))
                inner = match.group(2)
                if repeat_count < 1:
                    return False
                return _validate(inner)

            exposures = [s.strip() for s in seq_str.split(",")]
            for exp in exposures:
                parts = exp.split("/")
                if len(parts) != 3:
                    return False
                try:
                    n_repeat = int(parts[0])
                    filt = parts[1]
                    exp_time_str = parts[2].replace(",", ".")
                    if cmd == "FOCUS" and exp_time_str.lower() == "a":
                        pass
                    else:
                        exp_time = float(exp_time_str)
                        if exp_time < 0:
                            return False
                    if n_repeat < 1:
                        return False
                    if allowed_filters and filt not in allowed_filters:
                        return False
                except ValueError:
                    return False
            return True

        result["seq"] = _validate(seq_str)
        return result

    # ------------------------- static seq (BC forwarder) ---------------

    @staticmethod
    def validate_seq(obs: dict, allowed_filters=None) -> dict:
        """Backward-compatible static wrapper around the new seq module."""
        seq = obs.get("seq")
        if seq is None:
            return {}
        ok = _validate_seq(
            seq,
            command_name=obs.get("command_name"),
            allowed_filters=allowed_filters,
        )
        return {"seq": ok}


# ---------------------------------------------------------------------------
# Schema introspection helpers
# ---------------------------------------------------------------------------

def _collect_properties(schema: dict) -> dict:
    """Walk a schema and collect every declared property name + its sub-schema.

    Handles `properties`, nested `allOf`, `oneOf`, `anyOf`, and resolved `$ref`
    targets (only those already inlined). Used for type coercion and for the
    `allowed` hint returned by `validate_ob`.
    """
    props: dict = {}

    def visit(node: Any) -> None:
        if not isinstance(node, dict):
            return
        for k, v in node.get("properties", {}).items():
            props.setdefault(k, v)
        for keyword in ("allOf", "anyOf", "oneOf"):
            for sub in node.get(keyword, []):
                visit(sub)

    visit(schema)
    return props


def _collect_required(schema: dict, command_name: str | None) -> list[str]:
    """Collect `required` entries for the matching command branch, if any."""
    if not command_name:
        return sorted(schema.get("required", []))

    required: set[str] = set(schema.get("required", []))

    def visit(node: Any) -> None:
        if not isinstance(node, dict):
            return
        cn = node.get("properties", {}).get("command_name", {})
        if isinstance(cn, dict) and cn.get("const") == command_name:
            required.update(node.get("required", []))
        for keyword in ("allOf", "anyOf", "oneOf"):
            for sub in node.get(keyword, []):
                visit(sub)

    visit(schema)
    return sorted(required)


def _coerce_union(value: Any, prop: dict) -> Any:
    """Best-effort coercion for schemas that declare `oneOf: [int, str]` etc."""
    for alt in prop.get("oneOf", []) or prop.get("anyOf", []) or ():
        t = alt.get("type")
        if t == "integer":
            try:
                return int(value)
            except (ValueError, TypeError):
                continue
        if t == "number":
            try:
                return float(value)
            except (ValueError, TypeError):
                continue
        if t == "string":
            return str(value)
    return value
