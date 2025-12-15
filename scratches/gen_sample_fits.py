from astropy.io import fits
import numpy as np


imgsize = (128, 128)
# Empty data creation 128x128
data = np.zeros(imgsize, dtype=np.uint16)

# Empty header creation
header = fits.Header()

# Manual adding entries with comments
header['SIMPLE']   = (True, 'conforms to FITS standard')
header['BITPIX']   = (16, 'array data type')
header['NAXIS']    = (2, 'number of array dimensions')
header['NAXIS1']   = imgsize[0]
header['NAXIS2']   = imgsize[1]
header['OCASTD']   = ('1.1.2   ', 'Version of the OCA FITS headers standard')
header['OBSERVAT'] = ('OCA     ', 'Observatory name')
header['OBS-LAT'] = (-24.598056, '[deg] Observatory latitude')
header['OBS-LONG'] = (-70.196389, '[deg] Observatory longitude')
header['OBS-ELEV'] = (2817, '[m] Observatory elevation')
header['ORIGIN']   = ('CAMK PAN', 'Institution created this FITS file')
header['TELESCOP'] = ('dumy    ', 'Telescope name')
header['PI']       = ('Pietrzyn', 'Name of the Principal Investigator')
header['SCIPROG']  = ('ALICJA', 'Name of the science program')
header['OBS-ID']   = ('zvVgRjP0.004.05', 'BLOCK-ID.SEQ-ID.LOOP-ID (8.3.2) base64.dec.dec')
header['UOBI']     = ('96e9403b', 'Unique observing block id')
header['DATE-OBS'] = ('2024-06-01T23:27:46.764435', 'DateTime of observation start')
header['JD']       = (2460463.477624588, 'Julian date of observation start')
header['RA']       = (190.5875, 'Requested field RA (epoch=EQUINOX)')
header['DEC']      = (-0.553333, 'Requested field DEC (epoch=EQUINOX)')
header['EQUINOX']  = (2000.0, '[yr] Equinox of equatorial coordinates')
header['RA_OBJ']   = ('', 'Program object RA')
header['DEC_OBJ']  = ('', 'Program object DEC')
header['RA_TEL']   = (190.90214460702, 'Telescope mount RA (epoch=now)')
header['DEC_TEL']  = (-0.6876271096903661, 'Telescope mount DEC (epoch=now)')
header['ALT_TEL']  = (60.45281583118317, '[deg] Telescope mount ALT')
header['AZ_TEL']   = (38.67967882171911, '[deg] Telescope mount AZ')
header['AIRMASS']  = ('', 'Airmass')
header['OBSMODE']  = ('', 'Observation mode')
header['FOCUS']    = (15376, 'Focuser position')
header['ROTATOR']  = (53.619049072265625, '[deg] Rotator position')
header['OBSERVER'] = ('', 'Observers who acquired the data')
header['IMAGETYP'] = ('science ', 'Image type')
header['OBSTYPE']  = ('science ', 'Observation type')
header['OBJECT']   = ('SA104   ', 'Object name')
header['OBS-PROG'] = ('', 'Name of the science project')
header['NLOOPS']   = (2, 'Number of all exposures in this sequence')  # not for BESO
header['LOOP']     = (1, 'Number of exposure within this sequence')  # not for BESO
header['FILTER']   = ('B       ', 'Filter')
header['EXPTIME']  = (60.0, '[s] Executed exposure time')
# header['SENSOR'] = ('DW936_BV', 'Instrument name')
header['INSTRUME'] = ('DW936_BV', 'Instrument name') # 'SENSOR' versus INSTRUMENT
header['T-CAM']    = (-60.42599868774414, '[deg C] Temperature - current ccd/cmos')
header['T-CAMSET'] = (-60, '[deg C] Temperature - set ccd/cmos')
header['XBINNING'] = (1, 'Ccd binx')
header['YBINNING'] = (1, 'Ccd biny')
header['READ-MOD'] = (2, 'Readout mode')
header['GAIN-MOD'] = (2, 'Gain mode')
header['GAIN']     = ('', '[e-/ADU] Gain')
header['RON']      = ('', '[e-/read] Readout noise')
header['SUBRASTR'] = ('', 'Subraster size')
header['SCALE']    = (0.5, '[arcsec/pixel] Image scale')
header['SATURATE'] = ('', 'Data value at which saturation occurs')
# header['WCSAXES']  = (2, 'Number of coordinate axes')
# header['CRPIX1']   = (796.577613831, 'Pixel coordinate of reference point')
# header['CRPIX2']   = (1728.89355469, 'Pixel coordinate of reference point')
# header['PC1_1']    = (2.69285403195E-05, 'Coordinate transformation matrix element')
# header['PC1_2']    = (0.000137333902689, 'Coordinate transformation matrix element')
# header['PC2_1']    = (-0.000137383411552, 'Coordinate transformation matrix element')
# header['PC2_2']    = (2.68694419723E-05, 'Coordinate transformation matrix element')
# header['CDELT1']   = (1.0, '[deg] Coordinate increment at reference point')
# header['CDELT2']   = (1.0, '[deg] Coordinate increment at reference point')
# header['CUNIT1']   = ('deg     ', 'Units of coordinate increment and value')
# header['CUNIT2']   = ('deg     ', 'Units of coordinate increment and value')
# header['CTYPE1']   = ('RA---TAN', 'TAN (gnomonic) projection + SIP distortions')
# header['CTYPE2']   = ('DEC--TAN', 'TAN (gnomonic) projection + SIP distortions')
# header['CRVAL1']   = (190.684148716, '[deg] Coordinate value at reference point')
# header['CRVAL2']   = (-0.489958879133, '[deg] Coordinate value at reference point')
# header['LONPOLE']  = (180.0, '[deg] Native longitude of celestial pole')
# header['LATPOLE']  = (-0.489958879133, '[deg] Native latitude of celestial pole')
# header['MJDREF']   = (0.0, '[d] MJD of fiducial time')
# header['RADESYS']  = ('FK5     ', 'Equatorial coordinate system')
header['BSCALE']   = 1   # not in int
header['BZERO']    = 32768
header['TRACKING'] = (True, 'Mount tracking:  T=on F=off')
header['GUIDING']  = (False, 'Guiding status:  1=on 0=off')
header['DOME-AZ']  = (321.4, '[deg] Dome azimuth, 0 is North')
header['DOME-SHT'] = (2, 'Dome: 0,2=open(ing) 1,3=close(ing) 4&5=err')
header['MIRR-COV'] = (1, 'Mirror cover: 0&4=unkn 1=closed 2=move 3=open')
header['FAN-MIRR'] = (False, 'Ventilator mirror status:  T=on, F=off')
header['FAN-DOME'] = (True, 'Ventilator dome status:  T=on, F=off')
header['ROT-MECH'] = (221.1, '[deg] Rotator mechanical position')
header['LAMP-FLT'] = (True, 'Dome flat light: 0=unknown 1=off 2=busy 3=on')
header['P-WS']     = (1024.6, '[hPa] Pressure - weather station')
header['RHUM-WS']  = (21.7, '[%] Relative humidity - weather station')
header['T-WS']     = (12.9, '[deg C] Temperature - weather station')
header['WIND-DIR'] = (34.8, '[deg] Wind direction - weather station, 0 is N')
header['WIND-AVG'] = (12.1, '[m/s] Wind average speed - weather station')
header['WIND-GUS'] = (14.7, '[m/s] Wind gust - weather station')

header['COMMENT'] = ""

# Save
fname = 'example_ocm_raw_1_1_2.fits'
hdu = fits.PrimaryHDU(data=data, header=header)
hdu.writeto(fname, overwrite=True)

# Check file
from astropy.wcs import WCS

with fits.open(fname) as hdul:
    header = hdul[0].header
    wcs = WCS(header)
    print(f"Successfully loaded WCS from {fname}")
    print(f"WCS info: {wcs}")
