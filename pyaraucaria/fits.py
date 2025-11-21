import os
from collections import OrderedDict
from typing import List, Dict, Union, Optional

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
    hdul.writeto(file_name, overwrite=overwrite)
    if dtyp=='sideint16':
        with fits.open(file_name) as hdul2:
            hdul2[0].header['BZERO'] = 32768
            hdul2.writeto(file_name, overwrite=True)

# Let's Follow the FitS standard version 4, as defined in
# https://fits.gsfc.nasa.gov/fits_standard.html
# https://fits.gsfc.nasa.gov/standard40/fits_standard40aa-le.pdf
# https://heasarc.gsfc.nasa.gov/docs/fcg/common_dict.html may be also useful

# Mirk PROP:
# APERTURE FOC_LEN CREATOR (soft.) OBS_ID MOON_RA MOON_DEC or MOONANGL? DATAMAX DATAMIN
# PIXSIZE1= 3.800000E+00 / Pixel Size 1 (microns) PIXSIZE2=
# PIERSIDE OFFSET ?


def fits_header(
    oca_std="1.1.2",
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
    pi_name='',
    sci_prog='',
    obs_id='',
    uobi='',
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
    saturate='',
    pierside='',
    flat_era='',
    zero_era='',
    dark_era='',
    test='',
    baseline_clamp='',
    sensor_compensation='',
    sensor_port="",
    vertical_shift_speed='',
    rotator_mech_pos='',
    tracking='',
    guiding='',
    dome_az='',
    dome_sht='',
    mirror_cover='',
    fan_mirr='',
    fan_dome='',
    lamp_flat='',
    press_ws='',
    humidity_ws='',
    temp_ws='',
    wind_dir_ws='',
    wind_avg_ws='',
    wind_gust_ws='',
    ):

    _header = OrderedDict({
        "OCASTD": (oca_std, "OCA FITS HDU standard version"),
        "OBSERVAT": (obs, "Observatory name"),
        "OBS-LAT": (obs_lat, f"[deg] Observatory longitude"),
        "OBS-LONG": (obs_lon, f"[deg] Observatory latitude"),
        "OBS-ELEV": (obs_elev, f"[m] Observatory elevation"),
        "ORIGIN": (origin, "Institution created this FITS file"),
        "TELESCOP": (tel_id, 'Telescope name'),
        "PI": (pi_name, 'Name of the Principal Investigator'),
        "SCIPROG": (sci_prog, 'Name of the science program'),
        "OBS-ID": (obs_id, 'BLOCK-ID.SEQ-ID.LOOP-ID (8.3.2) base64.dec.dec'),
        "UOBI": (uobi, 'Unique observing block id'),
        "DATE-OBS": (utc_now, "DateTime of observation start"),
        "JD": (jd, "Julian date of observation start"),
        "RA": (req_ra, "Requested field RA"),
        "DEC": (req_dec, "Requested field DEC"),
        "EQUINOX": (epoch, '[yr] Equinox of equatorial coordinates'),
        "RA_OBJ": (ra_obj, "Program object RA"),
        "DEC_OBJ": (dec_obj, "Program object DEC"),
        "RA_TEL": (tel_ra, "Telescope mount RA"),
        "DEC_TEL": (tel_dec, "Telescope mount DEC"),
        "ALT_TEL": (tel_alt, "[deg] Telescope mount ALT"),
        "AZ_TEL": (tel_az, "[deg] Telescope mount AZ"),
        "AIRMASS": (airmass, 'Airmass'),
        "OBSMODE": (obs_mode, "Observation mode"),  # values exampl.: "TRACKING, GUIDING, JITTER or ELSE"
        "OBSERVER": (observer, 'Observers who acquired the data'),
        "IMAGETYP": (image_type, 'Image type'), # values exampl.: zero, flat, dark, science, focus
        "OBSTYPE": (obs_type, 'Observation type'),  # values exampl.: science, test, calib, art
        "OBJECT": (object, 'Object name'),
        # "OBS-PROG": (obs_prog, 'Name of the science project'),
        "NLOOPS": (n_loops, 'Number of all exposures in this sequence'),
        "LOOP": (loop, 'Number of exposure within this sequence'),
        "FILTER": (filter, 'Filter'),
        "EXPTIME": (exp_time, "[s] Executed exposure time"),
        "INSTRUME": (instrume_name, 'Instrument name'),  # full instrument name, like: 'Andor iKon-L DW936_BV'
        "T-CAM": (ccd_temp, '[deg C] Temperature - current ccd/cmos'),
        "T-CAMSET": (set_temp, '[deg C] Temperature - set ccd/cmos'),
        "XBINNING": (binx, 'Ccd binx'),
        "YBINNING": (biny, 'Ccd biny'),
        "READ-MOD": (read_mod, 'Readout mode'),
        "GAIN-MOD": (gain_mod, 'Gain mode'),
        "GAIN": (gain, '[e-/ADU] Gain'),
        "RON": (r_noise, '[e-/read] Readout noise'),
        "SUBRASTR": (subraster, 'Subraster size'),
        "SCALE": (scale, '[arcsec/pixel] Image scale'),
        "SATURATE": (saturate, 'Data value at which saturation occurs'),
        "FLAT_ERA": (flat_era, 'FLAT images era, increases on changes'),
        "ZERO_ERA": (zero_era, 'ZERO images era, increases on changes'),
        "DARK_ERA": (dark_era, 'DARK images era, increases on changes'),
        "CCD-BLCL": (baseline_clamp, 'Baseline clamp'),
        "CCD-SCMP": (sensor_compensation, 'Sensor compensation'),
        "CCD-PORT": (sensor_port, 'Sensor port (amplifier)'),
        "CCD-VSSP": (vertical_shift_speed, '[us/pixel] Vertical shift speed'),
        "FOCUS": (focus, "Focuser position"),
        "ROTATOR": (rotator_pos, "[deg] Rotator position"),
        "ROT-MECH": (rotator_mech_pos, '[deg] Rotator mechanical position'),
        "PIERSIDE": (pierside, 'EQ mount side of pier'),
        "TRACKING": (tracking, 'Mount tracking:  T=on F=off'),
        "GUIDING": (guiding, 'Guiding status:  1=on 0=off'),
        "DOME-AZ": (dome_az, '[deg] Dome azimuth, 0 is North'),
        "DOME-SHT": (dome_sht, 'Dome: 0,2=open(ing) 1,3=close(ing) 4&5=err'),
        "MIRR-COV": (mirror_cover, 'Mirror cover: 0&4=unkn 1=closed 2=move 3=open'),
        "FAN-MIRR": (fan_mirr, 'Ventilator mirror status:  T=on, F=off'),
        "FAN-DOME": (fan_dome, 'Ventilator dome status:  T=on, F=off'),
        "LAMP-FLT": (lamp_flat, 'Dome flat light: 0=unknown 1=off 2=busy 3=on'),
        "P-WS": (press_ws, '[hPa] Pressure - weather station'),
        "RHUM-WS": (humidity_ws, '[%] Relative humidity - weather station'),
        "T-WS": (temp_ws, '[deg C] Temperature - weather station'),
        "WIND-DIR": (wind_dir_ws, '[deg] Wind direction - weather station, 0 is N'),
        "WIND-AVG": (wind_avg_ws, '[m/s] Wind average speed - weather station'),
        "WIND-GUS": (wind_gust_ws, '[m/s] Wind gust - weather station'),
        "TEST": (test, 'Whether it is test frame, 0 or 1'),
        "COMMENT": (comment, 'Comment'),

    })
    # vertical_shift_speed: 134.49  # us/pixel   # header: CCD-VSSP
    # #          vertical_shift_speed: 156.89  # us/pixel
    # baseline_clamp: False  # header: CCD-BLCL
    # sensor_port: "all"  # which amplifier to use: all, bottom-left, bottom-right, top-left, top-right # header: CCD-PORT
    # sensor_compenstation: False  # header: CCD-SCMP

    return _header


def fits_stat(
        array: Union[np.ndarray, List], size: Optional[int] = None,
        min: bool = True, max: bool = True, mean: bool = True, median: bool = True,
        std: bool = True, rms: bool = True , sigma_quantile: bool = True) -> Dict:
    """
    Main fits statistics
    :param array: list like image data, example: [[2, 3, 4], [1, 1, 1], [2, 3, 1]] or numpy array
    :param size: size of sample taken random (if size=None will calculate whole array)
    :param min: array minimum, if True will estimate
    :param max: array maximum, if True will estimate
    :param mean: array mean, if True will estimate
    :param median: array median, if True will estimate
    :param std: array standard deviation, if True will estimate
    :param rms: array standard deviation, if True will estimate
    :param sigma_quantile: array sigma quantile, if True will estimate
    :return: Dict[str, float] with results
    """

    result = {}

    if isinstance(array, list):
        array = np.array(array)

    if size is not None:
        try:
            array = array_random_subset_2d(array, size=size)
        except ValueError:
            pass

    if min:
        result['min'] = float(np.amin(array))

    if max:
        result['max'] = float(np.amax(array))

    if mean:
        result['mean'] = float(np.mean(array))

    if median:
        result['median'] = float(np.median(array))

    if std:
        result['std'] = float(np.std(array))

    if rms:
        result['rms'] = float(np.std(array))

    if sigma_quantile:
        result['sigma_quantile'] = float(np.median(array) - np.quantile(array, 0.159))

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
