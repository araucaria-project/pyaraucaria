#!/usr/bin/env python3
"""
Handles Objects.database and TAB.ALL files

Useful functions
----------------
- `parse_object_database`
- `parse_tab_all`
- `map_objects_aliases`

Command Line Usage
------------------

Part of oca-pipe. Enjoy, Mikolaj
"""

import logging
import os
import re

__version__ = 1.4

#  The Locations of databases files, add, update if needed
objects_database_locations = [
    '/work/corvus/software/fits-warehouse/Objects.database',                 # main location
    os.path.join(os.path.split(__file__)[0], 'databases/Objects.database'),  # local backup (secondary choice)
]

tab_all_locations = [
    '/work/corvus/ONGOING/TAB.ALL',                                          # main location
    os.path.join(os.path.split(__file__)[0], 'databases/TAB.ALL'),           # local backup (secondary choice)
]


def parse_objects_database(file_path=None, skip_errors=True, radec_decimal=False):
    """
    Reads `Objects.database` file, returns objects dictionary and alias mapping.

    Parameters
    ----------
    file_path : str, optional
        Path to `Objects.database` file. If not provided, objects_database_locations will be used
    skip_errors : bool
        If true, on an error parser tries to skip line instead of throwing exception

    Returns
    -------
    (dict, dict)
        First dict is an object info, the keys are object names, each entry contains dict with 'ra', 'dec', 'aliases' as
            in `Objects.database`
        Second dict id a mapping, where keys are aliases and values are corresponding object names, those aliases are
            the ones from `Objects.database` AND, ADDITIONALLY THEIR CANONIZED VERSIONS (see `canonized_alias()`)
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
                            try:  # try convert to float
                                obj[m['key']] = float(m['val'])
                            except ValueError:
                                obj[m['key']] = m['val'].strip('\"')
                        elif i == 0: obj['ra'] = ra_to_decimal(t) if radec_decimal else t         # no key=value
                        elif i == 1: obj['dec'] = dec_to_decimal(t) if radec_decimal else t
                        elif i == 2: obj['per'] = float(t.split('?')[0])  # remove trailing ?
                        elif i == 3: obj['hjd0'] = float(t.split('?')[0])
                    obj['aliases'] = []
                    aliases[name] = name
                    objects[name] = obj
                else:  # aliases
                    obj['aliases'] += list(tokens)
                    aliases.update({al: name for al in tokens})
                    aliases.update({canonized_alias(al): name for al in tokens})
            except Exception as e:
                logging.error('on line %d of %s: %s', ln, file_path, str(e))
                if not skip_errors:
                    raise e
    return objects, aliases


def parse_tab_all(file_path=None, skip_errors=True, radec_decimal=False):
    """
    Reads `TAB.ALL` file. Returns objects and groups dictionary

    Parameters
    ----------
    file_path : str, optional
        Path to `Objects.database` file. If not provided, tab_all_locations will be used
    skip_errors : bool
        If true, on an error parser tries to skip line instead of throwing exception

    Returns
    -------
    (dict, dict)
        First dict contains object infos, the keys are object names, each entry contains dict with 'ra', 'dec', 'group'...
        Second dict is a group infos, where keys are group names and each entry contains dict with group parameters
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
                    if group_closed:  #  begin new group
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
                            key = m['key']
                            val = m['val']
                        else:
                            try:
                                key = actual_columns[i]
                                val = t
                            except IndexError:
                                logging.warning('on line %d of %s: ignoring value %s of unknown column', ln, file_path, t)
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
                            logging.error('on line %d of %s: %s (ra dec)=(%s %s) can not be converted into decimal '
                                          'repr and will be removed',
                                          ln, file_path, name, obj.get('ra'), obj.get('dec'))
                            try: del obj['ra']
                            except LookupError: pass
                            try: del obj['dec']
                            except LookupError: pass
                        except LookupError:
                            pass
                    objects[name] = obj
            except Exception as e:
                logging.error('on line %d of %s: %s', ln, file_path, str(e))
                if not skip_errors:
                    raise e

    return objects, groups


def map_objects_aliases(tab_all_objects: dict, aliases: dict):
    """
    Returns tab_all dict, but with keys exchanges to corresponding aliases keys, and subset (subdict?) of
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
        od_name = aliases.get(k)
        if od_name is None:
            orphans[k] = o
        else:
            mapped[od_name] = o
    return mapped, orphans


def lookup_objects(*objects, tab_all=None, objects_database=None, radec_decimal=False):
    """
    Returns all the information found in Objects.database and TAB.ALL for specified objects

    Object names are subject to alias mapping. For each objects the key 'ta' contains information from
    TAB.ALL (if found), and key 'od' from Objects.database (if found)

    Parameters
    ----------
    *objects
        One or more objects names to lookup
    objects_database : str, optional
        Path to custom TAB.ALL
    tab_all : str, optional
        Path to custom Objects.database
    radec_decimal : bool

    Returns
    -------
    dict of dict
        For each parameter returns all the information found in Objects.database and TAB.ALL
    """
    dbase = ObjetsDatabase(tab_all=tab_all, objects_database=objects_database, radec_decimal=radec_decimal)
    return dbase.lookup_objects(*objects)


class ObjetsDatabase(object):

    def __init__(self, tab_all=None, objects_database=None, skip_errors=True, radec_decimal=False):
        self.objects_database = objects_database
        self.tab_all = tab_all
        self.tab_all_objects, self.tab_all_groups = parse_tab_all(
            file_path=tab_all, skip_errors=skip_errors, radec_decimal=radec_decimal)
        self.objects_database_objects, self.objects_database_aliases = parse_objects_database(
            file_path=objects_database, skip_errors=skip_errors, radec_decimal=radec_decimal)
        self.tab_all_objects_mapped, orphans = map_objects_aliases(
            tab_all_objects=self.tab_all_objects, aliases=self.objects_database_aliases)
        self.tab_all_objects_mapped.update(orphans)

    def lookup_objects(self, *objects):
        """Returns dictionary with all available properties of aliases given as parameters

        Example
        -------
        >>> od = ObjetsDatabase()
        >>> od.lookup_objects('lmc169_5:84583', 'SMC09')
        {'lmc169_5:84583': {'name': 'LMC37', 'ra': '05:29:48.11', 'dec': '-69:35:32.1', 'aliases': ['LMC-T2CEP-136', 'pole3', 'lmc169_5_84583']}, 'SMC09': {'name': 'SMC09', 'ra': '00:43:37.1', 'dec': '-73:26:25.4', 'aliases': ['smc_sc3-63371']}}
        """
        ret = {}
        for o in objects:
            info = {}
            a = self.resolve_alias(o)
            try:
                info.update(self.objects_database_objects[a])
            except LookupError:
                pass
            try:
                info['taball'] = self.tab_all_objects_mapped[a]
            except LookupError:
                pass
            ret[o] = info
        return ret

    def resolve_alias(self, alias):
        """Returns corresponding object ID or `alias` itself if there is no mapping"""
        try:
            return self.objects_database_aliases[canonized_alias(alias)]
        except LookupError:
            return alias

    def get_object_properties_aliases(self, object, include_canonized=False):
        """For a given object-id (not alias, resolve first), returns tuple `(properties, aliases)`

        `properties` is a depth 1 (flat) dict with all properties extracted from `objects.database` and `TAB.ALL`,
        (`TAB.ALL` overwrites `objects.database`)

        `aliases` is a list of aliasses as defined in `objects.database` plus canonized versions if `include_canonized`

        Parameters
        ----------
        object : str
            ID of the object (not alias!)
        include_canonized : bool
            If true caononized versions of aliases are added to returned `aliases`

        Example
        -------
        >>> od = ObjetsDatabase(radec_decimal=True)
        >>> od.get_object_properties_aliases('CEP03', include_canonized=False)
        ({'name': 'LMC-CEP-2532', 'ra': 84.018667, 'dec': -70.032111, 'per': 2.03534862, 'hjd0': 4507.8, 'group': 'lmcpuls2', 'mI': 15.74, 'mV': 17.3, 'pa': '3.73829,-1.19662,0.47108,-0.00065,0.11096', 'pb': '9.49448,-2.87812,1.03983,-0.55972,0.35256'}, ['LMC-CEP-2532', 'Cep-2532', 'Cep_2532', 'LMC176.8-48147', 'LMC176.8_48147', 'lmc_cep2532', 'ceph2532', 'ceph2532_18', 'ceph2532_19', 'ceph2532_28', 'ceph2532_29', 'ceph2532_30', 'ceph2532_31', 'lmc-cep-2532_25', 'lmc-cep-2532_26', 'lmc-cep-2532_29', 'cep_2532_01', 'cep_2532_02'])

        """
        info = self.objects_database_objects[object]
        try:
            aliases = info.pop('aliases')
            if include_canonized:
                can = {canonized_alias(a) for a in aliases} - set(aliases)
                aliases += list(can)
        except LookupError:  # no aliases
            aliases = []
        try:
            taball = self.tab_all_objects_mapped[object]
            info.update(taball)
        except LookupError:
            pass

        return info, aliases



    def lookup_group(self, group, include_members=True):
        """Lookup for TAB.ALL defined group of objects"""
        grp = self.tab_all_groups[group]
        if include_members:
            objects = [o for o in self.tab_all_objects_mapped if o.get('group', None) == group]
            grp['objects'] = self.lookup_objects(*objects)
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




_transl_table = str.maketrans({'_': '-', ' ': '-', '.': '-', ':': '-'})


def canonized_alias(alias: str):
    """"""
    return alias.translate(_transl_table).lower()


_sexadec_parser = re.compile(r'(?P<sign>[+\-])?(?P<A>\d\d?)[ :\-hH](?P<B>\d\d?)[ :\-mM](?P<C>\d\d?(?:\.\d*)?)')


def _parse_sexadec(sexadec : str):
    v = _sexadec_parser.match(sexadec)
    if v is None:
        raise ValueError(f'{sexadec} can not be converted to decimal representation')
    if v['sign'] == '-':
        sign = -1.0
    else:
        sign = 1.0
    return sign, float(v['A']), float(v['B']), float(v['C'])


def ra_to_decimal(hms):
    try:
        val = float(hms)
    except ValueError:
        sign, h, m, s = _parse_sexadec(hms)
        val = sign * ((((s / 60) + m) / 60) + h) / 24 * 360
    return round(val, 6)


def dec_to_decimal(dms):
    try:
        val = float(dms)
    except ValueError:
        sign, d, m, s = _parse_sexadec(dms)
        val = sign * ((((s / 60) + m) / 60) + d)
    return round(val, 6)


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
    ret = lookup_objects(*args.object, tab_all=args.tab_all, objects_database=args.objects_database,
                         radec_decimal=args.decimal)

    if args.yaml:
        import yaml
        try:
            print(yaml.safe_dump(ret,sort_keys=False, default_flow_style=None, width=500))
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

