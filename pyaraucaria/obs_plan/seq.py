"""Semantic validation for the `seq=` mini-language.

Grammar (informal):

    seq     := group (',' group)*
    group   := count '/' filter '/' exposure
             | count 'x' '(' seq ')'
    count   := <positive int> | 'INF'
    filter  := <identifier, whitelist-checked by caller>
    exposure := <non-negative float>
             |  'a'     (auto-exposure; only allowed for FOCUS and SKYFLAT)

The outer `pattern` in `types/seq.yaml` only checks gross structure. This
function is the semantic layer: filter whitelist, exposure-time sanity,
repeat-count positivity, and context-sensitive `a` handling.
"""
from __future__ import annotations

import re
from collections.abc import Iterable

_MULTIPLIER_RE = re.compile(r"^(\d+|INF)x\((.+)\)$")
_AUTO_EXP_COMMANDS = frozenset({"FOCUS", "SKYFLAT"})


def validate_seq(
    seq_str: str,
    *,
    command_name: str | None = None,
    allowed_filters: Iterable[str] | None = None,
) -> bool:
    """Return True if `seq_str` is a semantically valid exposure sequence.

    Parameters
    ----------
    seq_str:
        The raw value of the `seq=` kwarg.
    command_name:
        Used to decide whether `a` (auto-exposure) is allowed.
    allowed_filters:
        Optional whitelist of filter names. If given, every filter in the
        sequence must be a member.
    """
    if not isinstance(seq_str, str):
        return False

    allowed = set(allowed_filters) if allowed_filters else None
    auto_ok = command_name in _AUTO_EXP_COMMANDS

    return _validate(seq_str.strip(), allowed, auto_ok)


def _validate(seq: str, allowed: set[str] | None, auto_ok: bool) -> bool:
    m = _MULTIPLIER_RE.match(seq)
    if m:
        count_tok, inner = m.group(1), m.group(2)
        if count_tok != "INF" and int(count_tok) < 1:
            return False
        return _validate(inner, allowed, auto_ok)

    for exp in (s.strip() for s in seq.split(",")):
        parts = exp.split("/")
        if len(parts) != 3:
            return False
        count_s, filt, exp_s = parts
        try:
            n = int(count_s)
        except ValueError:
            return False
        if n < 1:
            return False
        if allowed is not None and filt not in allowed:
            return False
        if exp_s.lower() == "a":
            if not auto_ok:
                return False
        else:
            try:
                t = float(exp_s.replace(",", "."))
            except ValueError:
                return False
            if t < 0:
                return False
    return True
