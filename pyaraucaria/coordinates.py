"""
Handles fast transformation of decimal (float) and sexagesimal (str) coordinates

Compatible with python 2.7 and 3.6+

Useful functions
----------------
- `ra_to_decimal`
- `dec_to_decimal`
- `ra_to_sexagesimal`
- `dec_to_sexagesimal`

Licence: MIT
"""

import re


_sexagesimal_parser = re.compile(r'(?P<sign>[+\-])?(?P<A>\d\d?)[ :\-hH](?P<B>\d\d?)[ :\-mM](?P<C>\d\d?(?:\.\d*)?)')


def parse_sexagesimal(sexagesimal):
    # type: (str) -> (int, int, int, float)
    """Returns quadruple: (sign, dec/h, min, sec)"""
    v = _sexagesimal_parser.match(sexagesimal)
    if v is None:
        raise ValueError('{} can not be converted to decimal representation'.format(sexagesimal))
    v = v.groupdict()
    if v['sign'] == '-':
        sign = -1
    else:
        sign = 1
    return sign, int(v['A']), int(v['B']), float(v['C'])


def format_sexagesimal(deg, multiplier, sign, sep=':'):
    # type: (float, float, bool, str) -> str
    if len(sep) == 1:
        sep = sep * 2

    sig = 1 if deg > 0 else -1
    h = sig * multiplier * deg
    m = (h * 60) % 60
    s = (m * 60) % 60
    return '{}{:02d}{}{:02d}{}{:06.3f}{}'.format(
        ('+' if sig > 0 else '-') if sign else '',
        int(h),
        sep[0],
        int(m),
        sep[1],
        s,
        sep[2] if len(sep) == 3 else ''
    )


def ra_to_decimal(hms):
    # type: (str) -> float
    try:
        val = float(hms)
    except ValueError:
        sign, h, m, s = parse_sexagesimal(hms)
        val = sign * ((((s / 60) + m) / 60) + h) / 24 * 360
    return round(val, 6)


def dec_to_decimal(dms):
    # type: (str) -> float
    try:
        val = float(dms)
    except ValueError:
        sign, d, m, s = parse_sexagesimal(dms)
        val = sign * ((((s / 60) + m) / 60) + d)
    return round(val, 6)


def ra_to_sexagesimal(deg, sep=':'):
    # type: (float, str) -> str
    return to_hourangle_sexagesimal(deg, sep=sep, sign=False)


def dec_to_sexagesimal(deg, sep=':'):
    # type: (float, str) -> str
    return to_degminsec_sexagesimal(deg, sep=sep, sign=True)

def to_hourangle_sexagesimal(deg, sep=':', sign=False):
    # type: (float, str, bool) -> str
    return format_sexagesimal(deg, 24.0/360, sign=sign, sep=sep)

def to_degminsec_sexagesimal(deg, sep=':', sign=True):
    # type: (float, str, bool) -> str
    return format_sexagesimal(deg, 1, sign=sign, sep=sep)


