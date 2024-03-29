import os
from typing import List, Tuple, Dict
import numpy
from astropy.io import fits


class Focus:

    METHODS = ["rms", "rms_quad", ]

    @staticmethod
    def build_focus_sharpness_coef(list_a: List, deg: int, focus_keyword: str,
                                   crop: int, focus_list: List = None, range: List or None = None) -> Tuple:
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

            # mean = numpy.mean(data)
            # median = numpy.median(data)
            rms = numpy.std(data)
            sharpness = rms
            focus_list_ret.append(float(focus))
            sharpness_list_ret.append(sharpness)
            hdu.close()

        coef = numpy.polyfit(x=focus_list_ret, y=sharpness_list_ret, deg=deg)

        return focus_list_ret, sharpness_list_ret, coef


    @staticmethod
    def calculate(fits_path: str, focus_keyword: str = "FOCUS", focus_list: List = None,
                  crop: int = 10, method: str = "rms_quad", range: List or None = None) -> Tuple[int, Dict]:
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

            focus_list_ret, sharpness_list_ret, coef = Focus.build_focus_sharpness_coef(list_a=list_a,
                                                                                        deg=deg,
                                                                                        focus_keyword=focus_keyword,
                                                                                        crop=crop,
                                                                                        focus_list=focus_list,
                                                                                        range=range)

        # ##### RMS with parabolic fit ######
        elif method == "rms":
            deg = 2
            if len(list_a) < deg:
                raise ValueError(f"for {method} method at least 2 focus positions are required")

            focus_list_ret, sharpness_list_ret, coef = Focus.build_focus_sharpness_coef(list_a=list_a,
                                                                                        deg=deg,
                                                                                        focus_keyword=focus_keyword,
                                                                                        crop=crop,
                                                                                        focus_list=focus_list,
                                                                                        range=range)

        else:
            focus_list_ret = []
            sharpness_list_ret = []
            coef = numpy.array([])

        a = numpy.max(focus_list_ret)
        b = numpy.min(focus_list_ret)
        x = numpy.linspace(start=a, stop=b, num=1000)
        y = numpy.polyval(coef, x)
        k = numpy.argmax(y)
        max_sharpness_focus = x[k]

        if numpy.abs(numpy.max(sharpness_list_ret) - numpy.min(sharpness_list_ret)) < 5:
            status = "to small sharpness range"
        elif y[0] >= y[k] or y[-1] >= y[k]:
            status = "wrong range"
        else:
            status = "ok"

        calc_metadata = {"status": status, "poly_coef": coef, "focus_values": focus_list_ret,
                         "sharpness_values": sharpness_list_ret}

        return max_sharpness_focus, calc_metadata
