import os
from typing import List, Tuple, Dict, Optional
import numpy
import numpy as np
from astropy.io import fits
import scipy
from scipy.optimize import curve_fit

from pyaraucaria.ffs import FFS


class Focus:

    METHODS = ["rms", "rms_quad", "fwhm", "laplacian", "lorentzian"]

    @staticmethod
    def fwhm(
            array, saturation: float, threshold: float = 20,
            kernel_size: float = 6, fwhm: float = 4) -> float | None:
        ffs = FFS(image=array)
        coo, _ = ffs.find_stars(threshold=threshold, kernel_size=kernel_size, fwhm=fwhm)
        if len(coo) > 3:
            fwhm_x, fwhm_y = ffs.fwhm(saturation=saturation)
            if fwhm_x and fwhm_y:
                return np.sqrt(fwhm_x ** 2 + fwhm_y ** 2)
            else:
                return None
        else:
            return None

    @staticmethod
    def lorentzian(x, a, x0, gamma, z):
        return a * (gamma ** 2) / ((x - x0) ** 2 + gamma ** 2) + z

    @staticmethod
    def build_focus_sharpness_coef(
            list_a: List, focus_keyword: str,
            crop: int, focus_list: List = None, range: Optional[List] = None, measurement: str = "rms") -> Tuple:

        focus_list_ret = []
        sharpness_list_ret = []

        for my_iter, f_file in enumerate(list_a):
            hdu = fits.open(f_file)
            hdr = hdu[0].header
            if focus_list is None:
                focus = hdr[focus_keyword]

            else:
                focus = focus_list[my_iter]
            if range is not None:
                if focus < range[0] or focus > range[1]:
                    continue
            data = hdu[0].data

            edge_rows = int(data.shape[0] * float(crop) / 100.)
            edge_cols = int(data.shape[1] * float(crop) / 100.)

            data = data[edge_rows:-edge_rows, edge_cols:-edge_cols]

            if measurement == 'rms':
                sharpness = numpy.std(data)
            elif measurement == 'laplacian':
                try:
                    laplacian = scipy.ndimage.laplace(data)
                    sharpness = float(np.var(laplacian))
                except ValueError:
                    sharpness = None
            elif measurement == 'fwhm':
                fwhm =  Focus.fwhm(
                    array=data,
                    saturation = data.max() * 0.8,
                )
                if fwhm is not None:
                    sharpness = 50 - fwhm
                else:
                    sharpness = None
            else:

                raise ValueError
            if sharpness:
                focus_list_ret.append(float(focus))
                sharpness_list_ret.append(float(sharpness))
            hdu.close()

        return focus_list_ret, sharpness_list_ret


    @staticmethod
    def calculate(fits_path: str, focus_keyword: str = "FOCUS", focus_list: List = None,
                  crop: int = 10, method: str = "rms_quad", range: Optional[List] = None) -> Tuple[int, Dict]:
        """
        Function to calculate the focus position of maximum sharpness for a given FITS files.

        Parameters:
        fits_path (str): The path to the FITS files directory or a list with FITS files.
        focus_keyword (str, optional): FIST file header keyword to retrive focus encoder position. Default: "FOCUS".
        focus_list (list or None, optional): A list of focus values to use for the calculation. If None, the focus
        values will be extracted from the FITS header. Defaults to None.
        crop (int, optional): The amount of pixels to crop from the edges of each image. Defaults to 10.
        method (str, optional): The method to use for calculating sharpness. Can be "rms", "rms_quad".
         Defaults to "rms_quad".

        Returns:
        tuple: (ax_sharpness_focus, calc_metadata).
        * ax_sharpness_focus: focus encoder value for maximum sharpness
        * calc_metadata: Dictionary of metadata has the following keys:
        - poly_coef: A NumPy array of coefficients for the polynomial fit used to calculate sharpness.
        - focus_values: A list of focus values used for the calculation.
        - sharpness_values: A list of corresponding sharpness values for each focus value.
        """

        if method not in Focus.METHODS:
            raise ValueError(f"Invalid method {method}")
        if isinstance(fits_path, list):
            if not all(os.path.isfile(file) for file in fits_path):
                raise ValueError(f"Invalid list with fits {fits_path}")
            list_a = fits_path
        else:
            if os.path.isdir(fits_path):
                list_a = [os.path.join(fits_path, f) for f in os.listdir(fits_path) if ".fits" in f]
            else:
                raise ValueError(f"{fits_path} is not valid dir")

        if focus_list:
            if not isinstance(fits_path, list):
                raise TypeError(f"if focus_list is provided, FITS files must be provided as a list")
            elif len(list_a) != len(focus_list):
                raise ValueError(f"focus_list and fits_files_list must have same length")

        # ##### RMS with quadratic fit ######
        if method == "rms_quad":
            deg = 4
            if len(list_a) < deg:
                raise ValueError(f"for {method} method at least 4 focus positions are required")

            focus_list_ret, sharpness_list_ret = Focus.build_focus_sharpness_coef(
                list_a=list_a,
                focus_keyword=focus_keyword,
                crop=crop,
                focus_list=focus_list,
                range=range,
                measurement='rms'
            )

            coef = numpy.polyfit(x=focus_list_ret, y=sharpness_list_ret, deg=deg)

            # ##### RMS with quadratic fit ######
        elif method == "lorentzian":
            deg = 4
            if len(list_a) < deg:
                raise ValueError(f"for {method} method at least 4 focus positions are required")

            focus_list_ret, sharpness_list_ret = Focus.build_focus_sharpness_coef(
                list_a=list_a,
                focus_keyword=focus_keyword,
                crop=crop,
                focus_list=focus_list,
                range=range,
                measurement='rms'
            )

            p0 = [
                max(sharpness_list_ret) - min(sharpness_list_ret),
                numpy.median(focus_list_ret),
                10,
                min(sharpness_list_ret)
            ]
            coef, _ = curve_fit(Focus.lorentzian, focus_list_ret, sharpness_list_ret, p0=p0)


        # ##### RMS with parabolic fit ######
        elif method == "rms":
            deg = 2
            if len(list_a) < deg:
                raise ValueError(f"for {method} method at least 2 focus positions are required")

            focus_list_ret, sharpness_list_ret = Focus.build_focus_sharpness_coef(
                list_a=list_a,
                focus_keyword=focus_keyword,
                crop=crop,
                focus_list=focus_list,
                range=range,
                measurement='rms'
            )
            coef = numpy.polyfit(x=focus_list_ret, y=sharpness_list_ret, deg=deg)

        # ##### FWHM with quadratic fit ######
        elif method == "fwhm":
            deg = 4
            if len(list_a) < deg:
                raise ValueError(f"for {method} method at least 4 focus positions are required")

            focus_list_ret, sharpness_list_ret = Focus.build_focus_sharpness_coef(
                list_a=list_a,
                focus_keyword=focus_keyword,
                crop=crop,
                focus_list=focus_list,
                range=range,
                measurement='fwhm'
            )

            coef = numpy.polyfit(x=focus_list_ret, y=sharpness_list_ret, deg=deg)

        # ##### laplacian variation with quadratic fit ######
        elif method == "laplacian":
            deg = 4
            if len(list_a) < deg:
                raise ValueError(f"for {method} method at least 4 focus positions are required")

            focus_list_ret, sharpness_list_ret = Focus.build_focus_sharpness_coef(
                list_a=list_a,
                focus_keyword=focus_keyword,
                crop=crop,
                focus_list=focus_list,
                range=range,
                measurement="laplacian"
            )

            coef = numpy.polyfit(x=focus_list_ret, y=sharpness_list_ret, deg=deg)

        else:
            focus_list_ret = []
            sharpness_list_ret = []
            coef = numpy.array([])

        a = numpy.max(focus_list_ret)
        b = numpy.min(focus_list_ret)
        x = numpy.linspace(start=a, stop=b, num=200)
        if method == "lorentzian":
            y = Focus.lorentzian(x, *coef)
            max_sharpness_focus = int(np.round(coef[1]))
        else:
            y = numpy.polyval(coef, x)
            k = numpy.argmax(y)
            max_sharpness_focus = int(x[k])

        if numpy.abs(numpy.max(sharpness_list_ret) - numpy.min(sharpness_list_ret)) < 5:
            status = "to small sharpness range"
        elif max_sharpness_focus < min(focus_list_ret) or max_sharpness_focus > max(focus_list_ret):
            status = "wrong range"
        # elif y[0] >= y[k] or y[-1] >= y[k]:
        #     status = "wrong range"
        else:
            status = "ok"

        calc_metadata = {"status": status, "coef": coef, "focus_values": focus_list_ret,
                         "sharpness_values": sharpness_list_ret, "fit_x": x, "fit_y": y}

        return max_sharpness_focus, calc_metadata
