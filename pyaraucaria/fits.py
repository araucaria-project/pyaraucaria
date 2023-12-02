import os
from collections import OrderedDict
from astropy.io import fits
import numpy as np


def save_fits_from_array(array,
                         folder: str,
                         file_name: str,
                         header,
                         overwrite: bool = False,
                         dtyp: str = 'int32',
                         do_not_scale_image_data: bool = False,
                         ignore_blank: bool = False,
                         uint: bool = False,
                         scale_back: bool or None = None):
    """
    Save fits file from array to selected location.
    Parameters
    ----------
    array - list like image data, example: [[2, 3, 4], [1, 1, 1], [2, 3, 1]]
    folder - folder, where fits will be saved, example: '/home/fits/'
    file_name - file name without '.fits', example: 'file_no_2233'
    header - fits header in dict format, example: {"FITS_STD": "beta_1", "TEL": "iris"}
    overwrite - overwrite existing file, default=False
    dtyp - array type [str], example: 'int16', 'int32'
    do_not_scale_image_data - astropy PrimaryHDU  parameter
    ignore_blank - astropy PrimaryHDU  parameter
    uint - astropy PrimaryHDU  parameter
    scale_back - astropy PrimaryHDU  parameter
    """
    if os.path.splitext(file_name)[1] == "":
        file_name = f"{file_name}.fits"
    file_name = os.path.join(folder, file_name)

    hdr = fits.Header()

    if isinstance(header, dict):
        for n in header.keys():
            try:
                hdr[n] = header[n][0]
                hdr.comments[n] = header[n][1]
            except (LookupError, TypeError):
                hdr[n] = header[n]
    elif isinstance(header, fits.Header):
        hdr = header
    else:
        hdr["OCASTD"] = "No fits header provided"

    if dtyp=='int32':
        narray = np.array(array, dtype=np.int32)
    elif dtyp=='int16':
        narray = np.array(array, dtype=np.int16)
    elif dtyp=='sideint16':
        s_array = np.array(array) - 32768
        narray = np.array(s_array, dtype=np.int16)
    elif dtyp=='none':
        narray = array
    else:
        narray = array

    hdu = fits.PrimaryHDU(data=narray,
                          header=hdr,
                          do_not_scale_image_data=do_not_scale_image_data,
                          ignore_blank=ignore_blank,
                          uint=uint,
                          scale_back=scale_back)
    hdul = fits.HDUList([hdu])
    if dtyp == 'sideint16':
        hdul[0].header['BZERO'] = 32768
    hdul.writeto(file_name, overwrite=overwrite)

# Let's Follow the FitS standard version 4, as defined in
# https://fits.gsfc.nasa.gov/fits_standard.html
# https://fits.gsfc.nasa.gov/standard40/fits_standard40aa-le.pdf
# https://heasarc.gsfc.nasa.gov/docs/fcg/common_dict.html may be also useful

# Mirk PROP:
# APERTURE FOC_LEN CREATOR (soft.) OBS_ID MOON_RA MOON_DEC or MOONANGL? DATAMAX DATAMIN
# PIXSIZE1= 3.800000E+00 / Pixel Size 1 (microns) PIXSIZE2=
# PIERSIDE OFFSET ?


def fits_header(oca_std="BETA3",
                obs="OCA",
                obs_lat='',
                obs_lon='',
                obs_elev='',
                origin='CAMK PAN',
                tel_id='',
                utc_now='',
                jd='',
                req_ra='',
                req_dec='',
                epoch='',
                ra_obj='',
                dec_obj='',
                tel_ra='',
                tel_dec='',
                tel_alt='',
                tel_az='',
                airmass='',
                obs_mode='',
                focus='',
                rotator_pos='',
                observer='',
                image_type='',
                obs_type='',
                object='',
                obs_prog='',
                n_loops='',
                loop='',
                filter='',
                exp_time='',
                instrume_name='',
                ccd_temp='',
                set_temp='',
                binx='',
                biny='',
                read_mod='',
                gain_mod='',
                gain='',
                r_noise='',
                subraster='',
                comment='',
                scale='',
                saturate=''
                ):

    _header = OrderedDict({
        "OCASTD": (oca_std, "OCA FITS HDU standard version"),
        "OBSERVAT": (obs, "Observatory name"),
        "OBS-LAT": (obs_lat, f"[deg] Observatory longitude"),
        "OBS-LONG": (obs_lon, f"[deg] Observatory latitude"),
        "OBS-ELEV": (obs_elev, f"[m] Observatory elevation"),
        "ORIGIN": (origin, "Institution created this FITS file"),
        "TELESCOP": (tel_id, 'Telescope name'),
        "DATE-OBS": (utc_now, "DateTime of observation start"),
        "JD": (jd, "Julian date of observation start"),
        "RA": (req_ra, "Requested field RA"),
        "DEC": (req_dec, "Requested field DEC"),
        "EQUINOX": (epoch, "Requested RA DEC epoch"),
        "RA_OBJ": (ra_obj, "Program object RA"),
        "DEC_OBJ": (dec_obj, "Program object DEC"),
        "RA_TEL": (tel_ra, "Telescope mount RA"),
        "DEC_TEL": (tel_dec, "Telescope mount DEC"),
        "ALT_TEL": (tel_alt, "[deg] Telescope mount ALT"),
        "AZ_TEL": (tel_az, "[deg] Telescope mount AZ"),
        "AIRMASS": (airmass, 'Airmass'),
        "OBSMODE": (obs_mode, "Observation mode"),  # values exampl.: "TRACKING, GUIDING, JITTER or ELSE"
        "FOCUS": (focus, "Focuser position"),
        "ROTATOR": (rotator_pos, "[deg] Rotator position"),
        "OBSERVER": (observer, 'Observers who acquired the data'),
        "IMAGETYP": (image_type, 'Image type'), # values exampl.: zero, flat, dark, science, focus
        "OBSTYPE": (obs_type, 'Observation type'),  # values exampl.: science, test, calib, art
        "OBJECT": (object, 'Object name'),
        "OBS-PROG": (obs_prog, 'Name of the science project'),
        "NLOOPS": (n_loops, 'Number of all exposures in this sequence'),
        "LOOP": (loop, 'Number of exposure within this sequence'),
        "FILTER": (filter, 'Filter'),
        "EXPTIME": (exp_time, "[s] Executed exposure time"),
        # "DETSIZE": (det_size, 'Detector size'),  # ?
        "INSTRUME": (instrume_name, 'Instrument name'),  # full instrument name, like: 'Andor iKon-L DW936_BV'
        "CCD-TEMP": (ccd_temp, 'Ccd actual temperature'),
        "SET-TEMP": (set_temp, 'Ccd set temperature'),
        "XBINNING": (binx, 'Ccd binx'),
        "YBINNING": (biny, 'Ccd biny'),
        "READ-MOD": (read_mod, 'Readout mode'),
        "GAIN-MOD": (gain_mod, 'Gain mode'),
        "GAIN": (gain, '[e-/ADU] Gain'),
        "RON": (r_noise, '[e-/read] Readout noise'),
        "SUBRASTR": (subraster, 'Subraster size'),
        # "CCDSEC": (ccd_sec, 'Ccd section'),
        "COMMENT": (comment, 'Comment'),
        "SCALE": (scale, '[arcsec/pixel] Image scale'),
        "SATURATE": (saturate, 'Data value at which saturation occurs'),

    })

    return _header


def fits_stat(array, size: int or None = None, min: bool = True, max: bool = True,
              mean: bool = True, median: bool = True, std: bool = True):
    """
    Main fits statistics
    :param array: list like image data, example: [[2, 3, 4], [1, 1, 1], [2, 3, 1]]
    :param size: size of sample taken random (if size=None will calculate whole array)
    :param min: array minimum, should calculate this stat, default Yes
    :param max: array maximum, should calculate this stat, default Yes
    :param mean: array mean, should calculate this stat, default Yes
    :param median: array median, should calculate this stat, default Yes
    :param std: array standard deviation, should calculate this stat, default Yes
    :return: Dict[str, float]
    """

    result = {}
    narray = np.array(array)
    if size is not None:
        narray = array_random_subset_2d(narray, size=size)

    if min:
        result['min'] = np.amin(narray)

    if max:
        result['max'] = np.amax(narray)

    if mean:
        result['mean'] = np.mean(narray)

    if median:
        result['median'] = np.median(narray)

    if std:
        result['std'] = np.std(narray)

    return result


def array_random_subset_2d(array, size: int, replace: bool = False):
    """
    Func randomly selecting points from array
    :param array: 2d numpy array
    :param size: size of subset
    :param replace: with or without replacement
    :return: array random subset
    """
    nrows, ncols = array.shape
    row_indices = np.random.choice(nrows, size=size, replace=replace)
    col_indices = np.random.choice(ncols, size=size, replace=replace)
    return array[row_indices, col_indices]
