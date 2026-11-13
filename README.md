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
pyaraucaria = { git = "https://github.com/araucaria-project/pyaraucaria.git"}
```
**Warning** If your project uses `poetry` the installed versions of poetry have to be `>=2.0.0` for git depenndency to `pyaraucaria`.

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
`pyaraucaria.ffs` provides fast frame statistics and quality estimates: background map, star detection, frame FWHM/ellipticity, satellite-trail and diffraction-spike detection, and pixel masks (saturation, cosmic rays, lines, bad columns, hot pixels).

```python
from pyaraucaria.ffs import FFS

ffs = FFS(image, saturation=45000, autobin=True)
ffs.calc_frame_fwhm()
ffs.stats["frame"]["fwhm"], ffs.stats["frame"]["used_stars"]
ffs.lines          # detected trails (astropy Table)
ffs.masks          # boolean masks, original image size
```

> **Binning.** By default (`bin=1, autobin=False`) the image is **not** binned. `autobin=True` (bin chosen so the
> longer side is ≤ 2500 px) or `bin=2` or more speeds up large frames (4096×4096 `calc_frame_fwhm()`: bin 1 ~12 s /
> 1.6 GB RAM, bin 2 ~9 s / 1.3 GB). Binning is used internally only (background, line and star detection). Coordinates, FWHM, masks and maps are always returned in the
> original image size, and pixel parameters are given in original pixels. Example: `examples/ffs.ipynb` (quick FWHM,
> full pipeline, calibration frames).

### Focus
`pyaraucaria.focus` implements various telescope focusing algorithms (RMS, FWHM, Lorentzian, Laplacian).

### Ephemeris
`pyaraucaria.ephemeris` offers calculations for moon illumination, object visibility, and other ephemeris-related data using `astropy`.

`pyaraucaria.ephemeris` exposes a clean, object-oriented API for use in your own Python projects.

#### Initialization

```python
from astropy.coordinates import EarthLocation
from astropy.time import Time
import astropy.units as u
from ocacal import Sun, Moon, Stars

# Define Observer Location
loc = EarthLocation(lat=-24.6*u.deg, lon=-70.2*u.deg, height=2400*u.m)
# Or use a named site
# loc = EarthLocation.of_site('paranal')

```

#### Calculating Events (Altitude Crossings)

Find precise times when the Sun reaches specific altitudes.

```python
sun = Sun(loc)

# Check for horizon (0), -6, and -12 degrees
# Returns a list of event dictionaries sorted by time
events = sun.get_events_by_altitude([0, -6, -12], start_time=Time.now())

for e in events:
    print(f"Altitude {e['target_alt']}° at {e['time_utc']} (Az: {e['az']:.1f}°)")

```

#### Batch Star Processing

Efficiently calculate data for multiple stars using vectorization.

```python
# Define your catalog
catalog = [
    {'id': 'Vega', 'ra': 279.23, 'dec': 38.78},
    {'id': 'Deneb', 'ra': 310.35, 'dec': 45.28}
]

stars = Stars(loc, catalog)

# Get positions for the next 5 hours in 1-hour steps
times = Time.now() + [0, 1, 2, 3, 4, 5] * u.hour
ephemeris = stars.get_ephemeris(times)

# Result is a dictionary keyed by Star ID
for star_id, data in ephemeris.items():
    print(f"--- {star_id} ---")
    for point in data:
        print(f"Time: {point['time_utc']} | Alt: {point['alt']:.2f}°")

```


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
