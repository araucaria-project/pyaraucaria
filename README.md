# OCA Observatory and Araucaria Project Common Routines Library

This library is intended to group various common routines and commandline tools used in  Araucaria Project and OCA observatory software.

Routines of this library should contain external dependencies as limited as possible and be compatible with Python 3.6+.

## Installation

#### Basic install:

```bash
    $ pip install git+https:https://github.com/araucaria-project/pyaraucaria.git
```
or, [more modern and safer](https://adamj.eu/tech/2020/02/25/use-python-m-pip-everywhere/)
```bash
    $ python -m pip install git+https:https://github.com/araucaria-project/pyaraucaria.git
```

#### Developer install
For those who want to contribute
```bash
    $ git clone https://github.com/araucaria-project/pyaraucaria.git
    $ cd pyaraucaria
    $ pip install -e ./
```

#### Usage in your project
Add the following line to your `requirements.txt` file:
```requirements.txt
git+https:https://github.com/araucaria-project/pyaraucaria.git
```

## Routines

### Lookup Objects
Lookup for objects/targets parameters using one of its aliases.
Uses `Objects.database` nad `TAB.ALL` files.

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

