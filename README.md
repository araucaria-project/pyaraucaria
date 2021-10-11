# Common Routines of OCA Observatory and Araucaria Project

This library is intended to group various common routines and commandline tools used in Araucaria Project and OCA observatory software.

Routines of this library should contain external dependencies as limited as possible and be compatible with Python 3.6+.

The compatibility with Python 2.7 is also valuable.

## Installation

#### Basic install:

```bash
    $ pip install git+https://github.com/araucaria-project/pyaraucaria.git
```
or, [more modern and safer](https://adamj.eu/tech/2020/02/25/use-python-m-pip-everywhere/)
```bash
    $ python -m pip install git+https://github.com/araucaria-project/pyaraucaria.git
```

#### Developer install
For those who want to contribute
```bash
    $ git clone https://github.com/araucaria-project/pyaraucaria.git
    $ cd pyaraucaria
    $ pip install -e ./
```

#### Usage in your project
Add the following line to your `requirements.txt` file (and/or to `install_requires` section of your `setup.py` if you use one):
```requirements.txt
git+https://github.com/araucaria-project/pyaraucaria.git
```

## Routines

### Coordinates
`coordinates.py` contains dependency-free fast routines to parse/format sexagesimal coordinates.

### Lookup Objects
Lookup for objects/targets parameters using one of its aliases.
Uses `Objects.database` nad `TAB.ALL` files.

Example

```python
>>> from pyaraucaria.lookup_objects import ObjectsDatabase
>>> od = ObjectsDatabase()
>>> od.lookup_object('lmc105_8_11987')
{'name': 'CEP25', 
 'ra': '05:18:12.8', 
 'dec': '-71:17:15.4', 
 'per': 3.4050955, 
 'hjd0': 2160.55457, 
 'aliases': ['LMC-T2CEP-085', 'lmc105.8_11987', 'lmc105_8_11987', 'cepii_lmc105_8_11987'], 
 'hname': 'LMC-T2CEP-085'
}
```

#### Command Line
```bash
$ lookup_objects -j hd167003
{"hd167003": {"name": "HIP37", "ra": "18:14:43.3", "dec": "-33:08:41.8", "aliases": ["hd_167003", "hd167003"]}}
```
See also
```bash
$ lookup_objects --help
```

#### Optional dependencies
For YAML output (`-y` option), some `yaml` python package should be installed.

### Star Object Data Library
Set of routines for parsing `TAB.ALL`-like files

Example
```python
>>> from pyaraucaria import libobject
>>> ol = libobject.ObjectList('TAB.ALL')
>>> ol.get_object('AL_Dor').data
{'I': None,
 'K': None,
 'V': 7.8,
 'aop': 1.8771,
 'band': 'V',
 'comment': None,
 'dec': '-60:36:14',
 'ecc': 0.1952,
 'file': None,
 'group': 'hipp',
 'hjd0': 2452764.100149,
 'k1': 57.477,
 'k2': 57.253,
 'lc': 'I',
 'obstype': None,
 'per': 14.90519957,
 'phext': None,
 'ra': '04:46:52.2',
 'status': None,
 'v0': 11.836}
```

