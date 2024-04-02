#!/usr/bin/env python3
"""
Handles Objects.database and TAB.ALL files

Compatible with python 2.7 and 3.6+

Useful functions
----------------
- `parse_object_database`
- `parse_tab_all`
- `map_objects_aliases`

Licence: MIT
"""

import logging
import os
import re
import string
from pyaraucaria.coordinates import ra_to_decimal, dec_to_decimal
from pyaraucaria.libobject import ObjectList, Object

logger = logging.getLogger('lookup_ob')

__version__ = 2.0

objects_database_locations = [
    '/work/corvus/software/fits-warehouse/Objects.database',                 # main location
    os.path.join(os.path.split(__file__)[0], 'databases/Objects.database'),  # local backup (secondary choice)
]

tab_all_locations = [
    '/work/corvus/ONGOING/TAB.ALL',                                          # main location
    os.path.join(os.path.split(__file__)[0], 'databases/TAB.ALL'),           # local backup (secondary choice)
]


def name_canonizator(name: str) -> str:
    """Name canonization by removing any non-alphanumeric characters and converting to lower case

    This should be used for astronomical objects name comparision, and for database lookup and indices.
    """
    # return re.sub(r'\W+', '', name).lower()
    return re.sub(r'[^\w]|_', '', name).lower()



def parse_objects_database(file_path=None, skip_errors=True, radec_decimal=False):
    """
    Reads `Objects.database` file, returns objects dictionary and alias mapping.

    Parameters
    ----------
    file_path : str, optional
        Path to `Objects.database` file. If not provided, objects_database_locations will be used
    skip_errors : bool
        If true, on an error parser tries to skip line instead of throwing exception
    radec_decimal : bool
        Convert coordinate values to decimal degrees

    Returns
    -------
    (dict, dict)
        First dict is an object info, the keys are object names, each entry contains dict with 'ra', 'dec', 'aliases' as
            in `Objects.database`
        Second dict id a mapping, where keys are aliases and values are corresponding object names, those aliases are
            the ones from `Objects.database`
    """
    if file_path is None:
        for fp in objects_database_locations:
            if os.path.exists(fp):
                file_path = fp
                break

    objects = {}
    aliases = {}
    with open(file_path) as fd:
        for ln, line in enumerate(fd):
            try:
                line = line.split('#')[0]  # remove trailing comments
                stripped = line.strip()
                if not stripped or stripped[0] == '#':
                    continue
                tokens = re.findall(r'\S+\s*=\s*\"[^\"]*\"|[^\s\"]+', stripped)  # like .split(' ') but handles quoted "
                if not line[0].isspace():
                    name = tokens[0]
                    obj = {'name': name}
                    for i, t in enumerate(tokens[1:]):
                        if t == '#':
                            pass
                        m = re.match(r'(?P<key>\w+)\s*=\s*(?P<val>\S?.*)', t)  # search for key=value
                        if m is not None:
                            m = m.groupdict()
                            try:  # try convert to float
                                obj[m['key']] = float(m['val'])
                            except ValueError:
                                obj[m['key']] = m['val'].strip('\"')
                        elif i == 0: obj['ra'] = ra_to_decimal(t) if radec_decimal else t         # no key=value
                        elif i == 1: obj['dec'] = dec_to_decimal(t) if radec_decimal else t
                        elif i == 2: obj['per'] = float(t.split('?')[0])  # remove trailing ?
                        elif i == 3: obj['hjd0'] = float(t.split('?')[0])
                    obj['aliases'] = []
                    # aliases[name] = name
                    objects[name] = obj
                else:  # aliases
                    obj['aliases'] += list(tokens)
                    if tokens:
                        obj['hname'] = tokens[0]
                    aliases.update({al: name for al in tokens})
                    # aliases.update({canonized_alias(al): name for al in tokens})  ## only original aliases
            except FutureWarning as e:
                logger.error('on line %d of %s: %s', ln, file_path, str(e))
                if not skip_errors:
                    raise e
    return objects, aliases


def canonized_keys(mapping):
    """
    Returns same dictionary, but with the keys canonized
    """
    return {canonized_alias(k): v for k, v in mapping.items()}


def parse_tab_all(file_path=None, skip_errors=True, radec_decimal=False):
    """
    Reads `TAB.ALL` file. Returns objects and groups dictionary

    Parameters
    ----------
    file_path : str, optional
        Path to `Objects.database` file. If not provided, tab_all_locations will be used
    skip_errors : bool
        PARAMETER IGNORED IN THE CURRENT VERSION
    radec_decimal : bool
        Convert coordinate values to decimal degrees


    Returns
    -------
    (dict, dict)
        First dict contains object infos, the keys are object names, each entry contains dict with 'ra', 'dec', 'group',
        Second dict is a group infos, where keys are group names and each entry contains dict with group parameters
    """
    if file_path is None:
        for fp in tab_all_locations:
            if os.path.exists(fp):
                file_path = fp
                break

    objects = {}
    groups = {}
    ol = ObjectList(file_path)
    for ob_id, val in ol.object_list.items():
        obj = {'hname': ob_id}
        for k, v in val.data.items():
            if v is None:
                continue
            if type(v) not in [float, int, bool, str, list]:
                v = str(v)
            obj[k] = v
        obj['type'] = val.type
        groups[obj['group']] = {'type': val.type}

        if radec_decimal:
            try:
                obj['ra'] = ra_to_decimal(obj['ra'])
                obj['dec'] = dec_to_decimal(obj['dec'])
            except ValueError:
                logger.warning('file: %s, coordinates of %s: (ra dec)=(%s %s) can not be converted into decimal '
                              'representation and will be removed',
                              file_path, ob_id, obj.get('ra'), obj.get('dec'))
                try:
                    del obj['ra']
                except LookupError:
                    pass
                try:
                    del obj['dec']
                except LookupError:
                    pass
            except LookupError:  # no ra dec
                pass

        objects[ob_id] = obj

    return objects, groups


def _parse_tab_all_original(file_path=None, skip_errors=True, radec_decimal=False):
    """ Old version, new version uses libobject.py
    """
    if file_path is None:
        for fp in tab_all_locations:
            if os.path.exists(fp):
                file_path = fp
                break
    objects = {}
    groups = {}
    default_columns = ['ra', 'dec', 'per', 'hjd0']
    actual_columns = default_columns
    with open(file_path) as fd:
        group = {}
        group_closed = True
        for ln, line in enumerate(fd):
            try:
                line = line.split('#')[0]  # remove trailing comments
                stripped = line.strip()
                if not stripped or stripped[0] == '#':
                    continue
                m = re.match(r'^\s*@\s*(?P<key>\w+)\s*=\s*(?P<val>\S?.*)$', stripped)  # group setting?
                if m is not None:
                    m = m.groupdict()
                    if group_closed:  # begin new group
                        group = {}
                        actual_columns = default_columns
                        group_closed = False
                    try:
                        group[m['key']] = float(m['val'])
                    except ValueError:
                        group[m['key']] = m['val']
                    if m['key'] == 'group':
                        groups[m['val']] = group
                    elif m['key'] == 'columns':
                        actual_columns = m['val'].split('|')
                else:
                    group_closed = True
                    tokens = re.findall(r'\S+\s*=\s*\"[^\"]*\"|[^\s\"]+', stripped)  # .split(' ') alike
                    name = tokens[0]
                    obj = {
                        'name': name,
                        'group': group['group'],
                    }
                    for i, t in enumerate(tokens[1:]):
                        if t == '#':
                            pass
                        m = re.match(r'(?P<key>\w+)\s*=\s*(?P<val>\S?.*)', t)  # search for key=value
                        if m is not None:
                            m = m.groupdict()
                            key = m['key']
                            val = m['val']
                        else:
                            try:
                                key = actual_columns[i]
                                val = t
                            except IndexError:
                                logger.warning('on line %d of %s: ignoring value %s of unknown column',
                                                ln, file_path, t)
                                continue
                        try:  # try convert to float
                            obj[key] = float(val)
                        except ValueError:
                            obj[key] = val.strip('\"')
                    if radec_decimal:
                        try:
                            obj['ra'] = ra_to_decimal(obj['ra'])
                            obj['dec'] = dec_to_decimal(obj['dec'])
                        except ValueError:
                            logger.error('on line %d of %s: %s (ra dec)=(%s %s) can not be converted into decimal '
                                          'repr and will be removed',
                                          ln, file_path, name, obj.get('ra'), obj.get('dec'))
                            try: del obj['ra']
                            except LookupError: pass
                            try: del obj['dec']
                            except LookupError: pass
                        except LookupError:
                            pass
                    objects[name] = obj
            except FutureWarning as e:
                logger.error('on line %d of %s: %s', ln, file_path, str(e))
                if not skip_errors:
                    raise e

    return objects, groups


def map_objects_aliases(tab_all_objects, aliases):
    # type: (dict, dict) -> (dict, dict)
    """
    Returns tab_all dict, but with keys exchanged to corresponding aliases keys, and subset (subdict?) of
    tab_all_objects for which there is no entry in aliases

    Look ups for `tab_all.keys()` in `aliases`, returns mapped dict and not-found dict. Every enrtry from
    `tab_all_objects` go to one of returned dicts. If you want to map what can be mapped, and leave the rest,
    just combine returned dicts:
    ```
        tab_all_objects, _ = parse_tab_all()
        _, aliases         = parse_object_database()
        mapped, orphans    = map_objects_aliases(tab_all_objects, aliases)
        mapped.update(orphans)
    ```

    Parameters
    ----------
    tab_all_objects : dict
        Usually, returned by `parse_tab_all`
    aliases : dict
        Usually, returned by `parse_objects_database`

    Returns
    -------
    (dict, dict)
        The first dict is like `tab_all_objects` but with keys from aliases (see also `return_unknown`),
        The second is a subset tab_all_objects with objects not found in aliases
    """
    mapped = {}
    orphans = {}
    for k, o in tab_all_objects.items():
        od_name = aliases.get(canonized_alias(k))
        if od_name is None:
            o['name'] = k
            orphans[k] = o
        else:
            o['name'] = od_name
            mapped[od_name] = o
    return mapped, orphans


def lookup_object(name, **kwargs):
    """Returns dictionary with all available properties an object given by name or alias
    """
    tab_all = kwargs.get('tab_all', None)
    objects_database = kwargs.get('objects_database', None)
    radec_decimal = kwargs.get('objects_database', False)
    dbase = ObjectsDatabase(tab_all=tab_all, objects_database=objects_database, radec_decimal=radec_decimal)
    return dbase.lookup_object(name)

def lookup_objects(objects_list, **kwargs):
    """
    Returns all the information found in Objects.database and TAB.ALL for specified objects

    Object names are subject to alias mapping. For each objects the key 'ta' contains information from
    TAB.ALL (if found), and key 'od' from Objects.database (if found)

    Parameters
    ----------
    objects_list
        List of objects names to lookup
    objects_database : str, optional
        Path to custom TAB.ALL
    tab_all : str, optional
        Path to custom Objects.database
    radec_decimal : bool

    Returns
    -------
    dict of dict
        For each parameter returns all the information found in Objects.database and TAB.ALL

    Notes
    -----
    Each call for this method parses database files. Consider using the `ObjectsDatabase` object if you call
    it frequently
    """
    tab_all = kwargs.get('tab_all', None)
    objects_database = kwargs.get('objects_database', None)
    radec_decimal = kwargs.get('objects_database', False)
    dbase = ObjectsDatabase(tab_all=tab_all, objects_database=objects_database, radec_decimal=radec_decimal)
    return dbase.lookup_objects(objects_list)


class ObjectsDatabase(object):

    def __init__(self, tab_all=None, objects_database=None, skip_errors=True, radec_decimal=False):
        self.objects_database = objects_database
        self.tab_all = tab_all
        self.tab_all_objects, self.tab_all_groups = parse_tab_all(
            file_path=tab_all, skip_errors=skip_errors, radec_decimal=radec_decimal)
        self.objects_database_objects, self.objects_database_aliases = parse_objects_database(
            file_path=objects_database, skip_errors=skip_errors, radec_decimal=radec_decimal)
        self.all_canonized_aliases = canonized_keys(self.objects_database_aliases)
        self.all_canonized_aliases.update({canonized_alias(k): k for k in self.objects_database_objects.keys()})
        self.tab_all_objects_mapped, orphans = map_objects_aliases(
            tab_all_objects=self.tab_all_objects, aliases=self.all_canonized_aliases)
        self.tab_all_objects_mapped.update(orphans)
        self.all_canonized_aliases.update({canonized_alias(k): k for k in orphans.keys()})
        self.all_objects = self.lookup_objects(
            self.objects_database_objects.keys() | self.tab_all_objects_mapped.keys())

    def lookup_object(self, obj):
        """Returns dictionary with all available properties an object given by name or alias
        """
        return self.lookup_objects([obj])[obj]

    def lookup_objects(self, objects_list):
        """Returns dictionary with all available properties of multiple aliases given as parameters

        The keys of the dictionary, are the names from the `*objects` argument, the value is the dictionary
        of the object properties

        Example
        -------
        >>> od = ObjectsDatabase(skip_errors=False)
        >>> od.lookup_objects(['lmc169_5:84583', 'SMC09'])
        {'lmc169_5:84583': {'name': 'LMC37', 'ra': '05:29:48.11', 'dec': '-69:35:32.1', 'aliases': ['LMC-T2CEP-136', 'pole3', 'lmc169_5_84583']}, 'SMC09': {'name': 'SMC09', 'ra': '00:43:37.1', 'dec': '-73:26:25.4', 'aliases': ['smc_sc3-63371']}}
        """
        ret = {}
        for o in objects_list:
            info = {}
            a = self.resolve_alias(o)
            try:
                info.update(self.objects_database_objects[a])
            except LookupError:
                pass
            try:
                info.update(self.tab_all_objects_mapped[a])
            except LookupError:
                pass
            if not info:
                info = None
            ret[o] = info
        return ret

    def get_object(self, alias):
        """Returns libobject.Object instance for specified alias"""
        d = self.lookup_object(alias)
        return Object(d['name'], data=d)

    def resolve_alias(self, alias):
        """Returns corresponding object ID or `alias` itself if there is no mapping"""
        try:
            return self.all_canonized_aliases[canonized_alias(alias)]
        except LookupError:
            return alias

    # @property
    # def all_objects(self):
    #     """Returns directory of all objects"""
    #     all_keys =
    #     ret = self.objects_database_objects.copy()
    #     ret.update(self.tab_all_objects_mapped)
    #     return ret

    # def get_object_properties_aliases(self, object):
    #     """For a given object-id (not alias, resolve first), returns triple `(properties, aliases, canonized)`
    #
    #     `properties` is a depth 1 (flat) dict with all properties extracted from `objects.database` and `TAB.ALL`,
    #     (`TAB.ALL` overwrites `objects.database`)
    #
    #     `aliases` is a list of aliasses as defined in `objects.database` plus canonized versions if `include_canonized`
    #
    #     `canonized` are `aliases` after canonization (lower case and with separators mapped to dash: `-`)
    #
    #     Parameters
    #     ----------
    #     object : str
    #         ID of the object (not alias!)
    #
    #     Example
    #     -------
    #     >>> od = ObjectsDatabase(radec_decimal=True)
    #     >>> od.get_object_properties_aliases('LMC60')
    #     ({'name': 'LMC60', 'ra': 81.97875, 'dec': -69.655389}, ['lmc_sc3-79892'], ['lmc-sc3-79892'])
    #
    #     """
    #     info = self.objects_database_objects[object]
    #     try:
    #         aliases = info.pop('aliases')
    #     except LookupError:  # no aliases
    #         aliases = []
    #     can = [canonized_alias(a) for a in aliases]
    #     try:
    #         taball = self.tab_all_objects_mapped[object]
    #         info.update(taball)
    #     except LookupError:
    #         pass
    #
    #     return info, aliases, can

    def lookup_group(self, group, include_members=True):
        """Lookup for TAB.ALL defined group of objects"""
        grp = self.tab_all_groups[group]
        if include_members:
            objects = [k for k, o in self.tab_all_objects_mapped.items() if o.get('group', None) == group]
            grp['objects'] = self.lookup_objects(objects)
        return grp

    _global_instances = {}

    @classmethod
    def get_instance(cls, tab_all=None, objects_database=None, skip_errors=True, radec_decimal=False):
        """Get existing instance if possible"""
        # print('instances: ', cls._global_instances)
        try:
            return cls._global_instances[(tab_all, objects_database, skip_errors, radec_decimal)]
        except LookupError:
            i = cls(tab_all=tab_all, objects_database=objects_database,
                    skip_errors=skip_errors, radec_decimal=radec_decimal)
            cls._global_instances[(tab_all, objects_database, skip_errors, radec_decimal)] = i
            return i

    def only_in_objects_database(self):
        return [k for k, o in self.objects_database_objects.items() if self.tab_all_objects_mapped.get(k)]


ObjetsDatabase = ObjectsDatabase  # Backward compatibility, there was a typo....

try:
    _transl_table = string.maketrans('_ .:', '----')  # python 2
except AttributeError:
    _transl_table = str.maketrans({'_': '-', ' ': '-', '.': '-', ':': '-'})  # python 3 only


def canonized_alias(alias):
    # type: (str) -> str
    """"""
    return alias.translate(_transl_table).lower()


def main():
    logging.basicConfig(level=logging.INFO, format='%(levelname)s %(asctime)s: %(message)s',
                        datefmt='%Y-%m-%dT%H:%M:%S %Z')

    """Command line entry"""
    import argparse

    parser = argparse.ArgumentParser(description='Lookups TAB.ALL and Objects.database for object(s). '
                                                 'The data from TAB.ALL has precedence. '
                                                 'Prints results in to parse JSON or YAML format'
                                                 '(pip install pyyaml for yaml)',
                                     epilog='Part of oca-pipe. Enjoy, Mikolaj'
                                     )
    parser.add_argument('object', nargs='+', type=str,
                        help='object name to look up')
    parser.add_argument('-y', '--yaml', action='store_true',
                        help=r'output complete info in YAML format (requires pyaml to be installed)')
    parser.add_argument('-j', '--json', action='store_true',
                        help=r'output complete info in JSON format')
    parser.add_argument('-c', '--coo', action='store_true',
                        help=r'output RA DEC coordinates also')
    parser.add_argument('-o', '--objects-database', type=str, metavar='FILE', default=None,
                        help=r'optional path to custom Objects.database')
    parser.add_argument('-t', '--tab-all', type=str, metavar='FILE', default=None,
                        help=r'optional path to custom TAB.ALL')
    parser.add_argument('-d', '--decimal',  action='store_true',
                        help=r'convert RA DEC to decimal representation')

    args = parser.parse_args()
    ret = lookup_objects(args.object, tab_all=args.tab_all, objects_database=args.objects_database,
                         radec_decimal=args.decimal)

    if args.yaml:
        import yaml
        try:
            print(yaml.safe_dump(ret, sort_keys=False, default_flow_style=None, width=500))
        except TypeError:
            print(yaml.safe_dump(ret))
    elif args.json:
        import json
        print(json.dumps(ret))
    elif args.coo:
        for q, d in ret.items():
            name = d.get('name', '(not found)')
            try:
                ra = d['ra']
                dec = d['dec']
            except LookupError:
                try:
                    ra = d['taball']['ra']
                    dec = d['taball']['dec']
                except LookupError:
                    ra = 'na'
                    dec = 'na'
            print(name, ra, dec)
    else:
        for q, d in ret.items():
            name = d.get('name', '(not found)')
            print(name)


if __name__ == "__main__":
    main()
