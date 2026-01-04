# Common Routines of OCA Observatory and Araucaria Project

`pyaraucaria` is a Python library providing a collection of common routines and command-line tools used in the Araucaria Project and OCA (Observatorio Cerro Armazones) observatory software. It aims to be a lightweight, dependency-focused toolkit for astronomical data processing and observatory operations.

The library requires Python 3.10+.

## Installation

### Basic Install

```bash
pip install git+https://github.com/araucaria-project/pyaraucaria.git
```

### Developer Install

```bash
git clone https://github.com/araucaria-project/pyaraucaria.git
cd pyaraucaria
uv sync --all-extras
```

### Usage in Your Project

Add to your `pyproject.toml` dependencies:
```toml
dependencies = [
    "pyaraucaria @ git+https://github.com/araucaria-project/pyaraucaria.git",
]
```
Or directly form PyPi, after checking versions (PyPi releases may lag behind GitHub):
```toml
dependencies = [
    "pyaraucaria>=2.11.0",
]
```

## Core Modules

### Coordinates
`pyaraucaria.coordinates` contains dependency-free, fast routines to parse and format sexagesimal coordinates (RA/Dec).

```python
from pyaraucaria.coordinates import ra_to_decimal, dec_to_sexagesimal
ra = ra_to_decimal('12:30:00')  # Returns 187.5
dec = dec_to_sexagesimal(-15.5)  # Returns '-15:30:00.000'
```

### Lookup Objects
Lookup for objects/targets parameters using one of its aliases.
Uses `Objects.database` and `TAB.ALL` files.

```python
from pyaraucaria.lookup_objects import ObjectsDatabase
od = ObjectsDatabase()
od.lookup_object('lmc105_8_11987')
```

### Date and Time
`pyaraucaria.date` handles conversion between Julian dates, datetime objects, and heliocentric corrections.

### FITS Operations
`pyaraucaria.fits` provides utilities for reading and writing FITS files, including header management and array saving.

### Fast FITS Statistics (FFS)
`pyaraucaria.ffs` provides optimized routines for star detection in images and basic image statistics (mean, median, noise estimation).

### Focus
`pyaraucaria.focus` implements various telescope focusing algorithms (RMS, FWHM, Lorentzian, Laplacian).

### Ephemeris
`pyaraucaria.ephemeris` offers calculations for moon illumination, object visibility, and other ephemeris-related data using `astropy` and `astroplan`.

### Airmass
`pyaraucaria.airmass` calculates airmass based on elevation using Kasten and Young's model.

### Reddening
`pyaraucaria.reddening` provides lookup for interstellar reddening based on coordinate databases (e.g., LMC).

### Dome Geometry
`pyaraucaria.dome_eq` calculates the required dome azimuth for telescopes on equatorial mounts.

### Observation Plan Parser
`pyaraucaria.obs_plan` contains a parser for custom observation plan formats using the `lark` grammar library.

## Command Line Tools

### `lookup_objects`
Query the object database from the command line.
```bash
lookup_objects -j hd167003
```

### `find_stars`
Perform star detection and calculate statistics on a FITS file.
```bash
find_stars path/to/image.fits gain=1.2 rn_noise=5.0
```

## Development and Testing

To run tests:
```bash
python -m unittest discover tests
```

## License
This project is licensed under the LGPL-3.0-or-later License.