import os

from astropy.io import fits
import numpy as np

def save_fits_from_array(array, folder: str, file_name: str, header, overwrite: bool = False):
    """
    Save fits file from array to selected location.
    Parameters
    ----------
    array - list like image data, example: [[2, 3, 4], [1, 1, 1], [2, 3, 1]]
    folder - folder, where fits will be saved, example: '/home/fits/'
    file_name - file name without '.fits', example: 'file_no_2233'
    header - fits header in dict format, example: {"FITS_STD": "beta_1", "TEL": "iris"}
    overwrite - overwrite existing file, default=False
    """
    if os.path.splitext(file_name)[1] == "":
        file_name = f"{file_name}.fits"
    file_name = os.path.join(folder, file_name)

    hdr = fits.Header()
    if isinstance(header, dict):
        for n in header.keys():
            hdr[n] = header[n]
    else:
        hdr["FITS_STD"] = "No fits header loaded"

    narray = np.array(array, dtype=np.int32)
    hdu = fits.PrimaryHDU(data=narray, header=hdr)
    hdul = fits.HDUList([hdu])
    hdul.writeto(file_name, overwrite=overwrite)
