# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Overview

`pyaraucaria` is a dependency-focused Python 3.10+ toolkit of common routines and CLI tools for the Araucaria Project and OCA (Observatorio Cerro Armazones) observatory software. It is a flat library of largely independent modules for astronomical data processing and observatory operations — there is no central framework tying the modules together; each addresses one domain (coordinates, dates, FITS, focus, ephemeris, object lookup, observation plans).

## Common Commands

Developer install uses `uv`:
```bash
uv sync --all-extras
```

Run the full test suite (tests use `unittest`, not pytest):
```bash
python -m unittest discover tests
```

Run a single test module / class / method:
```bash
python -m unittest tests.test_coordinate
python -m unittest tests.test_coordinate.TestCoordinate.test_ra_to_decimal
```

Console entry points (defined in `pyproject.toml` `[project.scripts]`):
```bash
lookup_objects -j hd167003          # pyaraucaria/lookup_objects.py:main
find_stars image.fits gain=1.2      # pyaraucaria/find_stars.py:main
```

## Architecture Notes

Things that require reading multiple files to understand:

- **Object lookup has a two-tier data-path fallback.** `lookup_objects.py` reads two flat-text databases (`Objects.database`, `TAB.ALL`), preferring shared observatory paths under `/work/corvus/...` and falling back to the bundled copies in `pyaraucaria/databases/`. `libobject.py` (legacy `ObjectList`/`Object`, authored separately, `from numpy import *` style) is the underlying parser/data model; `lookup_objects.py` wraps it and adds coordinate conversion. Name matching goes through `name_canonizator()` — strips non-alphanumerics and lowercases — so all name comparisons/indices are on the canonical form.

- **Observation-plan handling is a parse-then-validate pipeline.** `obs_plan/obs_plan_parser.py` (`ObsPlanParser`) parses the custom plan text format via an inline `lark` grammar (`LINE_GRAMMAR`) into command dicts; the Lark instance is built lazily and cached on the class. `ob_validator.py` (`ObsValidator`) then validates those dicts against JSON schemas, and can round-trip between flat "obdict" form and the command-string form (`convert_from_obdict` / `convert_to_obdict`). Schemas live in `pyaraucaria/schemas/*.yaml` and are loaded by name via `ObsValidator.load_schema()` (resolves relative to `__file__`). There are paired `base_*` and `tpg_*` schema/rules files. Grammar reference doc: `obs_plan/obsplan.md`.

- **Time/coordinate modules share a `_ensure_time` helper pattern.** Modules that accept flexible time inputs (e.g. `coordinates.py`, `ephemeris.py`) each define a local `_ensure_time` that normalizes `None`/`Time`/raw values into an astropy `Time`. `coordinates.py` is deliberately dependency-light and fast for sexagesimal↔decimal conversion; `ephemeris.py` is the astropy/astroplan-heavy layer (`Sun`, `Moon`, `Stars` classes, altitude-crossing events, vectorized batch ephemeris).

- **`ffs.py` (Fast FITS Statistics) is the largest, performance-sensitive module** — star detection and image statistics used by the `find_stars` CLI. See `examples/ffs_example.ipynb` and `examples/ffs_qc.ipynb`.

## Conventions

- Modules vary in age and style: newer modules use type hints and astropy; `libobject.py` is legacy (wildcard numpy import, module-level author/version headers). Match the style of the file you are editing rather than imposing one house style.
- Bundled package data (`pyaraucaria/databases/`, `pyaraucaria/schemas/`) is resolved relative to `__file__` and must be included in builds — it is listed under `[tool.hatch.build.targets.sdist]` in `pyproject.toml`. Add new data files there.
- Version lives in `pyproject.toml` (`version = ...`); bump it there for releases (PyPI may lag behind GitHub).
