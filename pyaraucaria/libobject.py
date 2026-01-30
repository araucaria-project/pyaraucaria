#################################################################
#                 Star Object Data Library                      #
author ="Bogumil Pilecki"                                       #
version="1.7"                                                  #
reldate="5 Aug 2021"                                           #
# Copyright (C) 2009-2014                                       #
# send comments and wishes to:                                  #
# E-mail: pilecki AT astrouw DOT edu DOT pl                     #
#################################################################
# Changes:
# version: 1.6
#   - Compatible with python 3 (Mikolaj, mkalusz@camk.edu.pl)
# version 1.7
#   - Prints changed to logging
#   - Lines commented-out does not finish group section
#################################################################
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.
#
#################################################################

import logging
from numpy import *
import sys, os, time, string
import ephem
from difflib import get_close_matches

log = logging.getLogger('libobject')

dpi = 2.*pi

def listize(sth,sep="+"):
    if type(sth) is str: sth=sth.split(sep)
    elif type(sth) is not list: sth=[sth]
    return sth

#def str2h(hms):
#    f = map(float, hms.split(":"))
#    val = f[0] + f[1]/60. + f[2]/3600.
#    return val

#def str2deg(dms):
#    if dms[0] == "-":
#        valsign = -1.0
#        dms = dms[1:]
#    f = map(float, dms.split(":"))
#    val = f[0] + f[1]/60. + f[2]/3600.
#    return val*valsign


#DAYS = {"month": 30., "M": 30., "week": 7., "W": 7., "2months": 60.}

#ALIASES = {"name": "id", "p": "per", "e": "ecc", "mag": "V"}
#alias2name = lambda name: ALIASES[name] if name in ALIASES else name


def get_jday_now():
    jdnow = ephem.julian_date(ephem.now())
    return jdnow

def get_jday_date(date):
    try:
        jd = ephem.julian_date(ephem.Date(date))
    except ValueError:
        log.critical("Wrong date format: %s", date)
        log.error(
            "Allowed is: yyyy/mm/dd hh:mm:ss with any number of trailing elements omitted (last element may be a floating point value, eg. 1980/05/30.75)")
        sys.exit(0)
    return jd

def get_ut_now():
    utnow = ephem.now()
    return utnow

def get_ut_date(date):
    try:
        ut = ephem.Date(date)
    except ValueError:
        log.critical("Wrong date format: %s", date)
        log.error(
            "Allowed is: yyyy/mm/dd hh:mm:ss with any number of trailing elements omitted (last element may be a floating point value, eg. 1980/05/30.75)")
        sys.exit(0)
    return float(ut)

def get_lt_now(force_tz=None):
    if force_tz is None:
        ltnow = ephem.localtime(ephem.now())
    else:
        ltnow = ephem.Date(float(get_ut_now())+force_tz/24.).datetime()
    return ltnow

def get_lt_date(date, force_tz=None):
    try:
        if force_tz is None:
            lt = ephem.localtime(ephem.Date(date))
        else:
            lt = ephem.Date(float(date)+force_tz/24.)
    except ValueError:
        log.critical("Wrong date format: %s", date)
        sys.exit(0)
    return lt

def ut2jday(date):
    jday = ephem.julian_date(date)
    return jday

def utvec2jday(date):
    jday = map(ephem.julian_date, date)
    return array(jday)



class ConfigurationSyntaxError(Exception):
    pass

"""
##############################################################
#############           Object List          #################
##############################################################
"""

class ObjectList():
    def __init__(self, fname=None, def_band = None, data_dir='data', group=None, sort=None):
        self.object_file = fname
        self.def_band = def_band
        self.data_dir = data_dir
        self.group = group
        self.def_sort = sort
        
        self.init_data()
        
        if fname is not None:
            self.load_objects(fname)
    
    def init_data(self):
        self.clear_defaults()
        self.clear_object_list()

    def get_objects(self, group=None, sort=None):
        if group is None:
            group = self.group
        
        if group is None:
            objects = self.object_list.values()
        else:
            objects = self.get_group(group)
        
        if sort is None:
            sort = self.def_sort
        if sort is not None:
            objects = self.sort_objects(objects, sort)
        
        return objects
    
    def load_objects(self, fname, mode="new"):

        self.clear_defaults()
        if mode == "new":
            self.clear_object_list()

        with open(fname, 'r') as fobj:
            for i, line in enumerate(fobj):
                if line.strip() == "":
                    self.defaults.pop("group", None) # remove 'group' keyword on empty line
                    continue
                if '#' in line:
                    line = line.split('#')[0]
                sline = line.strip()
                if len(sline) == 0:
                    continue
                elif sline[0] == '@':
                    self.parse_defaults(sline[1:], i+1)
                else:
                    # Comment: I'm not yet sure if all the objects should be read
                    # and only group is to be shown or we just read the group
                    # now: read the group only AND group selection when object list is requested
                    if self.group is not None:
                        if "group" not in self.defaults or self.group != self.defaults["group"]:
                            continue
                    oid, odata = self.parse_data(sline, i+1)
                
                    if oid in self.object_list.keys():
                        log.info("A duplicate ID detected: %s  - line %d ignored.", oid, (i + 1))
                        continue
                    what_to_load = []
                    if 'lc' in self.defaults:
                        what_to_load += ['lc']
                    if self.defaults.get('obstype',"sp")[0] == 'p':
                        what_to_load += ['phot']
                    else:
                        what_to_load += ['rv']     
                    obj = Object(oid, odata, self.defaults, load=what_to_load, data_dir=self.data_dir)
                    self.object_list[oid] = obj
                    self.object_order.append(obj)


    def clear_object_list(self):
        self.object_order = []
        self.object_list = {}

    
    def clear_defaults(self):
        self.defaults = {}

    def parse_defaults(self, defline, line_number):
        items = list(map(str.strip, defline.split('=')))
        if len(items) == 2:
            key, val = items
            if key == "sort" and self.def_sort is None:
                self.def_sort = val
                return
            if key == "columns":
                val = val.split('|')
            self.defaults[key] = val
        else:
            raise ConfigurationSyntaxError("error in defaults, line: %d"%line_number)
    
    def conv_spaces(self, line):
        conv = False
        cline = ""
        for i, c in enumerate(line):
            if c == '"':
                conv = not conv
            if conv and c == ' ':
                cline += '_'
            elif c != '"':
                cline += c
        return cline
    
    def parse_data(self, dataline, line_number):
        """ Convert string data line into a dictionary of keyword:value pairs. """

        try:
            if '"' in dataline:
                dataline = self.conv_spaces(dataline)

            items = dataline.split()

            assert len(items)>0, "object line should consist of at least ID field"

            # ID is necessary
            oid = items[0]
            items = items[1:]
            nitems = len(items)
            odata = {}

    #DEV        print "def:", self.defaults

            # data in predefined columns (only values there)
            if 'columns' in self.defaults:
                def_columns = self.defaults['columns']
                col_items = min([nitems, len(def_columns)])
                key_col = 0
                for i_col in range(col_items):
                    col_id = def_columns[i_col]
                    item = items[i_col]
                    if '=' not in item:
    #                    print col_id, item,
                        odata[col_id] = item
                        key_col += 1
                    else:
                        break
                keyitems = items[key_col:]
            else:
                keyitems = items

            # read data set by keyword=value pairs
            for keyitem in keyitems:
                keyval = keyitem.split('=')
                if len(keyval) == 2:
                    key, val = keyval
                    odata[key] = val
    #                print "*"+key, val,
                else:
                    raise ConfigurationSyntaxError("error in data, line: %d"%line_number)
        except Exception:
            log.critical("Error parsing data line %d: %s", line_number, dataline)
            raise

#        print

        return oid, odata

    def get_object(self, obj_id, human=False):
        # obj_id is really an Object
        if isinstance(obj_id, Object):
            return obj_id
        
        # get object IDs
        oids = self.object_list.keys()
        ovals = self.object_list.values()

        # obj_id is a string with part of or whole ID
        if type(obj_id) is str and not obj_id.isdigit():
            ids = map(str.lower, oids) #REMOVE [obj.data["id"].lower() for obj in self.objects]
            obj_id = obj_id.lower()
            match = get_close_matches(obj_id, ids, n=1, cutoff=0.3)
            if len(match)==0:
                log.warning("no ID similar to: %s found on a list", obj_id)
                return None
            obj_num = ids.index(match[0])
            return ovals[obj_num]

        # obj_id is a number of an object on a list
        if type(obj_id) is str: # and is a digit because it was already tested above
                obj_num = int(obj_id)
        elif type(obj_id) is int:
            obj_num = obj_id
        else:
            log.error("Bad ID specified.")
            return None

        if human: obj_num -= 1 
        
        nobj = len(self.object_list)
        if obj_num<0 or obj_num >= nobj:
            log.error("object number: %d not in range (1 - %d) for given list" % (obj_num + 1, nobj))
            return None
        return ovals[obj_num]

    def get_group(self, group):
        group_objects = []
        for obj in self.object_list.values():
            if group == obj.data['group']:
                group_objects += [obj]
        return group_objects
    
    def sort_objects(self, objects, psort):
        vals = []
        for obj in objects:
            if psort in obj.data:
                vals.append(obj.data[psort])
            else:
                vals.append(None)
        isort = argsort(array(vals))
        return [objects[i] for i in isort]
    
    
    def get_objdef_band(self):
        band_defs = {}
        for obj in self.object_list.values():
            band = obj.data.get('band', None)
            if band is not None:
                if band in band_defs.keys():
                    band_defs[band] += 1
                else:
                    band_defs[band] = 0
        if len(band_defs) > 0:
            i = argmax(band_defs.values())
            return band_defs.keys()[i]
        else:
            return None
        
    
    def get_best_band(self):
        mstats = {}
        for obj in self.object_list.values():
            mags = obj.get_all_mag()
            for m in mags.keys():
                if m in mstats.keys():
                    mstats[m] += 1
                else:
                    mstats[m] = 0
        if len(mstats) > 0:
            i = argmax(mstats.values())
            return mstats.keys()[i]
        else:
            return None
        
    
    def get_default_band(self, fallback='V'):
        if self.def_band is None:
            def_band = self.get_objdef_band()
            if def_band is None:
                def_band = self.get_best_band()
            if def_band is None:
                return fallback
            return def_band
        else:
            return self.def_band

    def set_group(self, group):
        self.group = group


"""
##############################################################
#############              Object            #################
##############################################################
"""

class Object():
    def __init__(self, oid=None, data={}, defaults={}, load=None, data_dir='data'):
        """ load - a list of keywords, tells what to load, for example: load=['rv','phot','lc']"""

        self._keywords = ["ra", "dec", "status", "comment", 'band', 'lc', 'group', 'obstype', 'file', 'phext']
        self._orbkeys = ["v0", "k1", "k2", "ecc", "aop"]
        self._typekeys = {"CONST": ['v0',"V","I","K"], "STD": ["V-I", "J-K"], "VAR": ["per", "hjd0"], "ECL": self._orbkeys, "PULS": ["pa", "pb", "pfi","sa","sb","sfi"]}
        self._aliases = {"P": "per", "vartype": "type", "remarks": "comment", 'stat': 'status', 'a0': 'v0'}
        self.bands = "UBVIRJHK"
        for band in self.bands:
            self._aliases["m"+band] = band

        # STD and VAR types are extened with magnitudes in different bands
        for vtype in ["STD","VAR"]:
            self._typekeys[vtype].extend(self._typekeys["CONST"])

        # ECL and PULS are extended with all VAR keys (including those from CONST)
        for vtype in ["ECL","PULS"]:
            self._typekeys[vtype].extend(self._typekeys["VAR"])

        # single keys that are to be splitted into many
        self._shortcuts = {"orb": self._orbkeys}

        # set of functions to initialize parameters
        self._fun_dict = {"hjd0": self.str2hjd}

        # convert value to list functions
        for key in ("pa", "pb", "pfi", "sa", "sb", "sfi"):
            self._fun_dict[key] = self._val2list
        # keep string values
        for key in ["id", "ra", "dec", "status", "comment", "band", 'lc', "group", 'obstype', 'phext', 'file']:
            self._fun_dict[key] = str

        if oid is not None:
            self.init_data(oid, data, defaults)

        self.rv_instr = {None: 0}
        self.rv = {"hjd": [], "ph": [], "rv1": []}
        self.phot_instr = {None: 0}
        self.phot = {"hjd": [], "ph": []}
        self.lcs = {}

        if 'data' in defaults:
            data_dir = defaults['data']

        # load all kind of data
        dloadfun = {'rv': self.load_rv, 'phot': self.load_phot, 'lc': self.load_lc}
        if oid is not None and load is not None:
            for key in load:
                dloadfun[key](data_dir=data_dir)


    def str2hjd(self, hjd):
        hjd = float(hjd)
        return self.fix_hjd(hjd)

    def fix_hjd(self, hjd):
        if hjd < 10000.:     hjd += 2450000.
        elif hjd < 100000.:  hjd += 2400000.
        return hjd
        
#    def _val2ra(self, ra):
#        if type(ra) is str:   ra = str2h(ra)
#        if ra>24:   ra /= 15.
#        return ra

#    def _val2dec(self, dec):
#        if type(dec) is str:   dec = str2deg(dec)
#        return dec

    def _val2list(self, val):
        if isinstance(val, str):
            return list(map(float, val.split(',')))
        else:
            return val

    def _str2val(self, key, sval):
        """ replaces strings as read from file into values used in the code (sometimes strings as well) """
        if sval is not None:
            if key in self._fun_dict:
                return self._fun_dict[key](sval)
            else:
                return float(sval)
        else:
            return None

    def _aliases2keys(self, idict, aliases):
        """ replace aliases with key names """
        for key in list(idict.keys()):
            if key in aliases:
                idict[aliases[key]] = idict[key]
    
    def _short2keys(self, idict, shortcuts):
        """ replace single shortcut keys with a set of keys """
        for key in shortcuts.keys():
            if key in idict:
                keys = shortcuts[key]
                vals = idict[key].split(',')
                # join the lists rejecting not used keys if len(vals) < len(keys)
                for k, v in zip(keys, vals):
                    idict[k] = v
            

    def init_data(self, oid, data={}, defaults={}):
        self.id = oid
        alldict = defaults.copy()
        alldict.update(data)
        self._aliases2keys(alldict, self._aliases)
        self._short2keys(alldict, self._shortcuts)
        self.type = alldict.get("type", "CONST")
        type_keys = self._keywords + self._typekeys[self.type]

        self.data = {}        
        for key in type_keys:
            self.data[key] = self._str2val(key, alldict.get(key))
        
#        print self.data
    
    def updata_data(self, data={}):
        type_keys = self._keywords + self._typekeys[self.type]
        for key in type_keys:
            self.data[key] = self._str2val(key, data.get(key))

    def get_fname(self, fname=None, data_dir="data", ext=[""]):
        if fname is not None:
            if os.path.isfile(fname): return fname
            log.error("Bad filename given: %s" % fname)
            return None
        for sid in [self.id, self.id.lower(), self.id.upper()]:
            for i_ext in ext:
                fname = data_dir + '/' + sid + '.' + i_ext
                if os.path.isfile(fname):
                    return fname
        return None

    def load_lc(self, data_dir="data"):
        if self.data['lc'] is None:
            return
        for b in self.data['lc']:
            fname = self.get_fname(data_dir=data_dir, ext=b)
            if fname is None: continue
            log.info("Loading:", fname)
            lc = {}
            lhjd, lmag, ler = [], [], []
            with open(fname) as f:
                for line in f:
                    sline = line.split()
                    obsvals = map(float, sline[:3])
                    lmag += [obsvals[1]]
    #                ler  += [obsvals[2]]        
                    lhjd += [self.fix_hjd(obsvals[0])]
            lc['hjd'] = array(lhjd)
            lc['ph'] = array(self.get_phase(lc['hjd']))
            lc['mag'] = array(lmag)
    #        self.lc["er"] = array(ler)
            self.lcs[b] = lc

    def get_phot_ext(self):
        phext = self.data['phext']
        if phext is not None:
            if phext[:1] == '.':
                return phext[1:]
            return phext
        return "dat"

    def load_phot(self, data_dir="data"):
        ext = self.get_phot_ext()
        if self.data['file'] is not None:
            fname = data_dir + '/' + self.data['file']
        else:
            fname = self.get_fname(None, data_dir,[ext])
            if fname is None: return
        if not os.path.isfile(fname):   return

        log.info("Loading:", fname)

        photdata = {'hjd':[], 'mag':[], 'err':[], 'mask':[], 'instr':[]}        

        with open(fname) as fphot:
            for n, line in enumerate(fphot):
                sline = line.split()
                if len(sline) == 0:
                    continue
                data = {'mag': 0.0, 'err': -1.0, 'mask': 0, 'instr': 0}
                data['hjd'] = self.str2hjd(sline.pop(0))

                instr = None        
                for sl, par in zip(sline, ['mag', 'err']):
                    if sl[0] in "0123456789-+.":
                        data[par] = float(sl)
                    else:
                        instr = sl
                        break

                if len(sline) == 3:
                    instr = sline[2]

                if instr is not None:
                    if instr not in self.phot_instr.keys():
                        self.phot_instr[instr] = len(self.phot_instr)
                    data['instr'] = self.phot_instr[instr]

                if data['hjd'] in photdata['hjd']:
                    log.warning("%s: observation date %s is not unique! Duplicate entry is ingored.", fname, data['hjd'])
                    continue
                for dkey in data:
                    photdata[dkey].append(data[dkey])
                
        for key in ['hjd', 'mag', 'err', 'mask']:
            self.phot[key] = array(photdata[key])

        self.phot["ph"] = array(self.get_phase(self.phot["hjd"]))
        self.phot["instr"] = array(photdata['instr'],int)
 

    def read_pulsrv_linelist(self, sline):
        """ Read and interpret a line of RV data file.
            mask = 0 - velocity not set (no velocity extracted from a measurement)
            mask = 1 - velocity set
            mask = 2 - velocity set, but inactive (object not seen in the data, or of poor quality)
            Returns:
                data, status
        """
    
        data = {'rv1': 0.0, 'err1': -1.0, 'mask1': 0, 'instr': 0}
        data['hjd'] = self.str2hjd(sline.pop(0))

        instr = None        
        for sl, par in zip(sline, ['rv1', 'err1']):
            if sl[0] in "0123456789-+.":
                data[par] = float(sl)
            else:
                instr = sl
                break

        if len(sline) == 3:
            instr = sline[2]

        if instr is not None:
            if instr not in self.rv_instr.keys():
                self.rv_instr[instr] = len(self.rv_instr)
            data['instr'] = self.rv_instr[instr]

        return data
            
        # OTHER VERSION IF THE ONE ABOVE FAILS        
        if len(sline) == 0:
            return data
        
        if sline[0][0] in "0123456789-+.":
            data['rv1'] = float(sline.pop(0))
        else:
            data['instr'] = sline.pop(0)

        if len(sline) == 0:
            return data
        
        if sline[0][0] in "0123456789-+.":
            data['err1'] = float(sline.pop(0))
        else:
            data['instr'] = sline.pop(0)
            
        if len(sline) > 0:
            data['instr'] = sline.pop(0)

        return data


    def read_eclrv_linelist(self, sline):
        """ Read and interpret a line of RV data file.
            mask = 0 - velocity not set (no velocity extracted from a measurement)
            mask = 1 - velocity set
            mask = 2 - velocity set, but inactive (object not seen in the data, or of poor quality)
            Returns:
                data, status
        """
    
        data = {'rv1': 0.0, 'rv2': 0.0, 'err1': -1.0, 'err2': -1.0, 'mask1': 1, 'mask2': 1, 'instr': 0}
        data['hjd'] = self.str2hjd(sline.pop(0))
        
        if len(sline) in [1,3,5]:
            instr = sline.pop()
            if instr not in self.rv_instr.keys():
                self.rv_instr[instr] = len(self.rv_instr)
            data['instr'] = self.rv_instr[instr]

        # if errors are provided: negative err value means inactive velocity measurement
        if len(sline) == 4:
            vals = map(float, sline)
            if vals[1] < 0.0:
                data['mask1']  = 2
            else:
                data['rv1']  = vals[0]
                data['err1'] = vals[1]
            if vals[3] < 0.0:
                data['mask2']  = 2
            else:
                data['rv2']  = vals[2]
                data['err2'] = vals[3]
        # if errors are NOT provided: '-' in velocity column means inactive velocity measurement
        elif len(sline) == 2:
            if sline[0] == '-':
                data['mask1'] = 2
            else:
                data['rv1']  = float(sline[0])
            if sline[1] == '-':
                data['mask2'] = 2
            else:
                data['rv2']  = float(sline[1])
        # velocities not set
        elif len(sline) == 0:
            data['mask1'] = 0
            data['mask2'] = 0
        else:
            return None

        return data

    def load_rv(self, fname=None, data_dir="data"):
        if self.data['file'] is not None:
            fname = data_dir + '/' + self.data['file']
        else:
            fname = self.get_fname(fname,data_dir,["rv", "dat"])
            if fname is None: return
        if not os.path.isfile(fname):   return

        log.info("Loading: %s", fname)

        if self.type == "ECL":
            rvdata = {'hjd':[], 'rv1':[], 'rv2':[], 'err1':[], 'err2':[], 'mask1':[], 'mask2':[], 'instr':[]}
            read_rv_list = self.read_eclrv_linelist
        elif self.type == "PULS":
            rvdata = {'hjd':[], 'rv1':[], 'rv2':None, 'err1':[], 'err2':None, 'mask1':[], 'mask2':None, 'instr':[]}
            read_rv_list = self.read_pulsrv_linelist         

        with open(fname) as frv:
            for n, line in enumerate(frv):
                sline = line.split()
                if len(sline) == 0:
                    continue
                data = read_rv_list(sline)
                if data is None:
                    log.error("Error reading data file %s in line %d" % (fname, n))
                    continue
                if data['hjd'] in rvdata['hjd']:
                    log.warning("%s: observation date is not unique! Duplicate entry is ingored.", fname, data['hjd'])
                    continue
                for dkey in data.keys():
                    rvdata[dkey].append(data[dkey])
                
        for key in ['hjd', 'rv1', 'rv2', 'err1', 'err2', 'mask1', 'mask2']:
            if rvdata[key] is not None:
                self.rv[key] = array(rvdata[key])

        self.rv["ph"] = array(self.get_phase(self.rv["hjd"]))
        self.rv["instr"] = array(rvdata['instr'],int)

#        print "DBG:", self.rv

    def set_star_ephem(self):
        ephem_data = "star,f|V|G2,%s,%s,9.99,2000"%(self.data["ra"],self.data["dec"])
        self.data["star"] = ephem.readdb(ephem_data)
        return self.data["star"]

    def get_phase(self, hjd=None, ut=None, hshift=0.0):
        if self.type in ["CONST", "STD"] or self.data['hjd0'] is None or self.data['per'] is None:
            return None
        if hjd is None:
            hjd = ephem.julian_date(ut)
        ph = (hjd+hshift/24. - self.data["hjd0"]) / self.data["per"]
        return ph%1.0

    def ph_is_good(self, ph):
        if self.type == "ECL":
            if ph>0.15 and ph<0.35 or ph>0.65 and ph<0.85: return True
        elif self.type == "PULS":
            if ph>0.8 or ph<0.1: return True
        return False

    def ph_is_vgood(self, ph):
        if self.type == "ECL":
            if ph>0.2 and ph<0.3 or ph>0.7 and ph<0.8: return True
        elif self.type == "PULS":
            if ph>0.9 or ph<0.05: return True
        return False

    def get_all_mag(self):
        mdata = {}
        for b in self.bands:
            if b in self.data.keys():
                if self.data[b] is not None:
                    mdata[b] = self.data[b]
    
    def get_mag(self, band=None):
        if band is None or band not in self.bands or band not in self.data:
            return None
        else:
            return self.data[band]

    # PULS RV CURVE
    def get_puls_rvcurve(self, ph):
        rvc = zeros(len(ph))
        if self.data['v0'] is not None:
            rvc += self.data['v0']
        if self.data['pfi'] is not None:
            for w, fun in [['pa',cos], ['pb',sin]]:
                if self.data[w] is not None:
                    if len(self.data['pfi']) == 1:
                        pfi = self.data['pfi']*len(self.data[w])
                    else:
                        pfi = self.data['pfi']
                    for i,(a,fi) in enumerate(zip(self.data[w], pfi)):
                        rvc += a*fun(dpi*(i+1)*(ph + fi))
        else:
            for w, fun in [['pa',cos], ['pb',sin]]:
                if self.data[w] is not None:
                    for i,val in enumerate(self.data[w]):
                        rvc += val*fun(dpi*(i+1)*ph)
        return rvc

    def get_puls_rvs(self, times, time_type="hjd"):
        if time_type == "ut":
                times = utvec2jday(times)

        ph = self.calc_phase(times) #'times' is hjd here
        rv = self.get_puls_rvcurve(ph)

        return ph, rv

    def get_phot4obs(self):
        if len(self.phot['hjd']) == 0:
            return empty(0)
        mag = self.phot["mag"].copy()
        mask = self.phot["mask"]
        ph, cmag = self.get_puls_rvs(times=self.phot["hjd"], time_type="hjd") # calculated photometry
        
        sel = (mask==0)
        mag[sel] = cmag[sel]
        
        return mag

    def get_puls_rv4obs(self):
        if len(self.rv['hjd']) == 0:
            return empty(0)
        rv1 = self.rv["rv1"].copy()
        m1 = self.rv["mask1"]
        ph, crv1 = self.get_puls_rvs(times=self.rv["hjd"], time_type="hjd") # calculated RV
        
        s1 = (m1==0)
        rv1[s1] = crv1[s1]
        
        return rv1

    # ECL RV CURVE
    def get_rvs(self, times, time_type="hjd"):
        if time_type == "ut":
                times = utvec2jday(times)

        ph = self.calc_phase(times) #'times' is hjd here
        if self.data.get("ecc") in [None, 0.0]:
            rv1 = self.data.get("v0", 0.0) - self.data.get("k1", 50.0)*sin(dpi*ph)
            rv2 = self.data.get("v0", 0.0) + self.data.get("k2", 50.0)*sin(dpi*ph)
        else:
            hjd0, per, warning = self.get_ephem()
            v0, k1, k2, ecc, aop = self.get_orbit()
            ptimes = times-per*((aop - pi/2.)/dpi)
            rv1 = star_rv(ptimes, v0, k1, ecc, aop, hjd0, per)
            rv2 = star_rv(ptimes, v0, -k2, ecc, aop, hjd0, per)
        return ph, rv1,rv2

    def get_rvs4obs(self):
        if len(self.rv['hjd']) == 0:
            return empty(0), empty(0)
        rv1 = self.rv["rv1"].copy()
        rv2 = self.rv["rv2"].copy()
        m1, m2 = self.rv["mask1"], self.rv["mask2"]
        ph, crv1,crv2 = self.get_rvs(times=self.rv["hjd"], time_type="hjd") # calculated RVs
        
        s1, s2 = (m1==0), (m2==0)
        rv1[s1] = crv1[s1]
        rv2[s2] = crv2[s2]
        
        return rv1, rv2

    def calc_phase(self, hjd):
        hjd0 = self.data['hjd0']
        per = self.data['per']
        if per is None:
            hjd0, per =  0.0, 1.0
        if hjd0 is None:
            hjd0 = 0.0
        phs = (hjd - hjd0) / per
        return phs

    def get_phases(self, start="now", extent=1.0, step=None, times=None, time_type="hjd"):
#        print start, extent, step
        if times is None:
            if type(start) is str:
                start = get_jday_now()
                if "yesterday".startswith(start):   start -= 1.0
                elif "tomorrow".startswith(start):  start += 1.0
            if step is None:
                if extent>7:
                    step =1.0
                else:
                    step = extent/24.
            times = arange(start, start + extent, step)
        else:
            if time_type == "ut":
                times = utvec2jday(times)
        phs = (times - self.data["hjd0"]) / self.data["per"]
        return phs % 1.0

    def get_ephem(self):
        if self.type in ["CONST", "STD"]:
            return 0.0, 1.0, self.type
        hjd0 = self.data['hjd0']
        per = self.data['per']
        if per is None:
            return 0.0, 1.0, "no period"
        if hjd0 is None:
            return 0.0, per, "no hjd0"
        return hjd0, per, None

    def get_data(self, key, none_def = None):
        if self.data[key] is None:
            return none_def
        else:
            return self.data[key]

    def get_orbit(self):
        v0 = self.get_data('v0', 0.0)
        k1 = self.get_data('k1', 50.)
        k2 = self.get_data('k2', 50.)
        ecc = self.get_data('ecc', 0.0)
        aop = self.get_data('aop', 0.0)
        return v0, k1, k2, ecc, aop
        


"""
################ OLD ############### OLD ################## OLD ####################
"""

class OldObject():

    def fix_phase_for_ecc(self, ph):
        if self.data["ecc"]>0.0:
            ph = get_ecc_anomaly(0.0, 0.0, self.data["ecc"], M=ph*dpi)/dpi
        return ph%1.0

    def fit_rvs(self, v0def=250, k1def=20, k2def=20):
        self.data["v0"] = v0def
        self.data["k1"] = k1def
        self.data["k2"] = k2def
        if len(self.rv["rv1"])>0:
            rv1cur = self.rv["rv1"]
            rv2cur = self.rv["rv2"]
            rv12 = r_[rv1cur,rv2cur]
            ymin = min(rv12)
            ymax = max(rv12[rv12<900])
            self.data["v0"] = (ymax+ymin)/2.
            if len(self.rv["rv1"])>1:
                self.data["k1"] = ymax - self.data["v0"]
                self.data["k2"] = self.data["v0"] - ymin

    def get_colors(self):
        mcol1 = self.rv["cn"]/10.
        mcol2 = self.rv["cn"]/10. + 0.5
        return mcol1, mcol2
        


    

###################################################################################        
###################################################################################


def star_rv(x, v0, k, e, aop, hjd0, P):
    ph=(x-hjd0)/P%1.0 #- (aop - pi/2.)/dpi   # phase 0.0 - phase of periastron = phase 0.0 at periastron
    v = get_true_anomaly(0.0, 0.0, e, M=ph*dpi)%dpi
    rv = (cos(v+aop) + e*cos(aop))/sqrt(1-e**2)
    return v0 + k*rv
    
def get_ecc_anomaly(P,t,e, tp=0.0, M=None):
    if M is None: M = get_mean_anomaly(P,t, tp)
    E=empty(len(M))
    for i,Mi in enumerate(M):
        E[i]=newton(lambda xE,yM,ye: xE-ye*sin(xE)-yM, Mi, lambda xE,yM,ye: 1.0-ye*cos(xE), args=(Mi,e))
    return E

def get_mean_anomaly(P,t,tp=0.0):
    return 2.0*pi*(t-tp)/P

def get_true_anomaly(P,t,e, tp=0.0, M=None):
    E=get_ecc_anomaly(P,t,e, tp,M)
    v = 2*arctan(sqrt((1+e)/(1-e))*tan(E/2))
    return v    
    
###############################           NUMERICAL METHODS             ##############################################
# Netwon-Raphson method taken from scipy.optimize.minpack.py
def newton(func, x0, fprime=None, args=(), tol=1.48e-8, maxiter=50):
    """Given a function of a single variable and a starting point,
    find a nearby zero using Newton-Raphson.

    fprime is the derivative of the function.  If not given, the
    Secant method is used.
    """
    if fprime is not None:
        p0 = x0
        for iter in range(maxiter):
            myargs = (p0,)+args
            fval = func(*myargs)
            fpval = fprime(*myargs)
            if fpval == 0:
                print("Warning: zero-derivative encountered.")
                return p0
            p = p0 - func(*myargs)/fprime(*myargs)
            if abs(p-p0) < tol:
                return p
            p0 = p
    else: # Secant method
        p0 = x0
        p1 = x0*(1+1e-4)
        q0 = func(*((p0,)+args))
        q1 = func(*((p1,)+args))
        for iter in range(maxiter):
            if q1 == q0:
                if p1 != p0:
                    print("Tolerance of %s reached" % (p1 - p0))
                return (p1+p0)/2.0
            else:
                p = p1 - q1*(p1-p0)/(q1-q0)
            if abs(p-p1) < tol:
                return p
            p0 = p1
            q0 = q1
            p1 = p
            q1 = func(*((p1,)+args))
    raise RuntimeError("Failed to converge after %d iterations, value is %s" % (maxiter,p))


if __name__ == "__main__":
    olist = ObjectList("EXAMPLE.DAT")
    print(len(olist.get_objects()), "objects loaded")




