import numpy as np
from scipy.signal import convolve2d, find_peaks
from scipy.ndimage.filters import maximum_filter
from scipy.ndimage import convolve

from astropy.stats import mad_std


#  pip install -e .

class FFS:
    """
        FFS (Fast Fits Statistics) library for star detection in an image, and basic statistics.

        Args:
            image (np.ndarray): The input image.
            gain (float, optional): Detector gain value. Used for noise estimation. Defaults to 1.
            rn_noise (float, optional): Detector readout noise. Used for noise estimation. Defaults to 0.

        Attributes:
            all Args
            min (float): The minimum value of the ADU.
            max (float): The maximum value of the ADU.
            mean (float): The mean value of the ADU.
            median (float): The median value of the ADU.
            rms (float): The root mean square value of the ADU.
            sigma_quantile (float): The sigma value calculated as the quantile 0.5 - 0.159 .
            noise (float): The noise calculated as the Poisson noise accounting gain and readout noise.

        Methods:
            find_stars(self,threshold=5.,method="sigma quantile",kernel_size=9,fwhm=2): Finds the stars in the image with specified noise calulation method.

              Args:
                  threshold (float, optional): The threshold value for star detection. Defaults to 5.
                  method (str, optional): The method used for determining the sigma value. Can be 'rms Poisson','rms','sigma quantile'. Defaults to "sigma quantile".
                  kernel_size (int, optional): The size of the Gaussian kernel. Defaults to 9.
                  fwhm (float, optional): FWHM value. Defaults to 2.

              Returns:
                  coo (np.ndarray): An sorted array of coordinates representing the positions of stars.
                  adu (np.ndarray): An sorted array of ADU values corresponding to the detected stars.

            fwhm(self,saturation=65000,radius=10,all_stars=False): Calculates the average fwhm for stars in the X and Y axis.

              Args:
                saturation (float): Saturation level above fwhm calculation will be ignored for a star. Defaults to 65000
                radius (int): Radius in which fwhm will be calculated. Defaults to 10
                all_stars (bool): If True, fwhm will be calculated for all stars.
                                  If False, only for 100 non saturated brightests. Defaults to False

              Returns:
                fwhm_x,fwhm_y (float,float): Median of fwhm for X and Y axis, respectively

              Attributes:
                fwhm_xarr (np.ndarray): array of fwhm in X axis for stars, ordered accordingly to star ADU
                fwhm_yarr (np.ndarray): array of fwhm in Y axis for stars, ordered accordingly to star ADU

        Example usage:
            stats = FFS(data,threshold=5,kernel_size=9,fwhm=6)
            sigma = stats.sigma_quantile
            p_noise = stats.noise
            coo,adu = stats.find_stars()
            fwhm_x,fwhm_u = stats.fwhm(saturation=50000,all_stars=True)
            fwhm_xarr = stats.fwhm_xarr
            fwhm_yarr = stats.fwhm_yarr
        """

    def __init__(self, image, gain=1., rn_noise=0.):
        self.image = np.transpose(image)
        self.gain = float(gain)
        self.rn_noise = float(rn_noise)
        self.stats = {}

    def mk_stats(self):
        self.min = np.min(self.image)
        self.max = np.max(self.image)
        self.mean = np.mean(self.image)
        self.median = np.median(self.image)
        self.rms = np.std(self.image)
        self.sigma_quantile = np.median(self.image) - np.quantile(self.image, 0.159)

        self.q_sigma_lower = np.median(self.image) - np.quantile(self.image, 0.159)
        self.q_sigma_upper = np.quantile(self.image, 0.841) - np.median(self.image)
        self.q_sigma = (np.quantile(self.image, 0.841) - np.quantile(self.image, 0.159)) / 2.

        self.noise = (self.median / self.gain + self.rn_noise) ** 0.5

        self.stats["min"] = self.min
        self.stats["max"] = self.max
        self.stats["mean"] = self.mean
        self.stats["median"] = self.median
        self.stats["rms"] = self.rms
        self.stats["q_sigma_lower"] = self.q_sigma_lower
        self.stats["q_sigma_upper"] = self.q_sigma_upper
        self.stats["q_sigma"] = self.q_sigma


    def find_stars(self, threshold=5, method="sigma quantile", kernel_size=30, fwhm=10):

        self.coo = []
        self.adu = []
        self.threshold = float(threshold)
        self.method = method
        self.fwhm_adopted = float(fwhm)
        self.kernel_size = int(kernel_size)
        self.kernel_sigma = float(fwhm) / 2.355
        self.kernel = self.gauss_kernel(self.kernel_size, self.kernel_sigma)

        if self.method == "rms Poisson":
            self.sigma = self.noise
        elif self.method == "rms":
            self.sigma = self.rms
        elif self.method == "sigma quantile":
            self.sigma = self.sigma_quantile
        else:
            raise ValueError(f"Invalid method type {self.method}")

        maska1 = self.image > self.median + self.threshold * self.sigma
        data2 = convolve2d(self.image, self.kernel, mode='same')
        maska2 = (data2 == maximum_filter(data2, 3))
        maska = np.logical_and(maska1, maska2)
        coo = np.argwhere(maska)
        self.stats["stars"] = {}
        if len(coo) > 1:
            self.coo = coo
            x, y = zip(*self.coo)
            val = self.image[x, y]
            sorted_i = np.argsort(val.astype(float))[::-1]
            sorted_coo = self.coo[sorted_i]
            sorted_val = val[sorted_i]
            self.coo = sorted_coo
            self.adu = sorted_val

            tmp = {"x": self.coo[:,0], "y": self.coo[:,1],"max_adu":self.adu}
            self.stats["stars"] = tmp

        return self.coo, self.adu

    def fwhm(self, saturation=65000, radius=10, all_stars=True):
        radius = int(radius)
        self.fwhm_xarr = []
        self.fwhm_yarr = []
        self.fwhm_x = None
        self.fwhm_y = None
        for i, tmp in enumerate(self.coo):
            if all_stars:
                i_max = len(self.adu)
            else:
                i_max = 100
            d1 = d2 = d3 = d4 = None
            if self.adu[i] < int(saturation) and i < i_max:
                x, y = self.coo[i]
                max_adu = self.adu[i]
                half_adu = (max_adu - self.median) / 2.

                if True:
                    line = self.image[x - radius:x + radius, y] - self.median - half_adu
                    line = self.image[x - radius + 1:x + 1, y] - self.median - half_adu
                    maska1, maska2 = line > 0, line < 0
                    pos, neg = line[maska1], line[maska2]
                    if len(pos) > 0 and len(neg) > 0:
                        lower, upper = max(neg), min(pos)
                        line = list(line)
                        lower_i, upper_i = line.index(lower), line.index(upper)
                        lower_adu, upper_adu = line[lower_i], line[upper_i]
                        d1 = radius - upper_i - np.abs(lower_adu) / (np.abs(lower_adu) + np.abs(upper_adu))

                    line = self.image[x:x + radius, y] - self.median - half_adu
                    maska1, maska2 = line > 0, line < 0
                    pos, neg = line[maska1], line[maska2]
                    if len(pos) > 0 and len(neg) > 0:
                        lower, upper = max(neg), min(pos)
                        line = list(line)
                        lower_i, upper_i = line.index(lower), line.index(upper)
                        lower_adu, upper_adu = line[lower_i], line[upper_i]
                        d2 = upper_i + 1 - np.abs(lower_adu) / (np.abs(lower_adu) + np.abs(upper_adu))

                    line = self.image[x, y - radius + 1:y + 1] - self.median - half_adu
                    maska1, maska2 = line > 0, line < 0
                    pos, neg = line[maska1], line[maska2]
                    if len(pos) > 0 and len(neg) > 0:
                        lower, upper = max(neg), min(pos)
                        line = list(line)
                        lower_i, upper_i = line.index(lower), line.index(upper)
                        lower_adu, upper_adu = line[lower_i], line[upper_i]
                        d3 = radius - upper_i - np.abs(lower_adu) / (np.abs(lower_adu) + np.abs(upper_adu))

                    line = self.image[x, y:y + radius] - self.median - half_adu
                    maska1, maska2 = line > 0, line < 0
                    pos, neg = line[maska1], line[maska2]
                    if len(pos) > 0 and len(neg) > 0:
                        lower, upper = max(neg), min(pos)
                        line = list(line)
                        lower_i, upper_i = line.index(lower), line.index(upper)
                        lower_adu, upper_adu = line[lower_i], line[upper_i]
                        d4 = upper_i + 1 - np.abs(lower_adu) / (np.abs(lower_adu) + np.abs(upper_adu))

            if d1 != None and d2 != None:
                dx = (d1 + d2)
            else:
                dx = 0

            if d3 != None and d4 != None:
                dy = (d3 + d4)
            else:
                dy = 0
            self.fwhm_xarr.append(dx)
            self.fwhm_yarr.append(dy)

        self.fwhm_xarr, self.fwhm_yarr = np.array(self.fwhm_xarr), np.array(self.fwhm_yarr)

        maska1 = self.fwhm_xarr == 0
        maska2 = self.fwhm_yarr == 0
        maska = maska1 & maska2
        fwhm_xarr = self.fwhm_xarr[~maska]
        fwhm_yarr = self.fwhm_yarr[~maska]

        if len(fwhm_xarr) > 2:
            self.fwhm_x = np.median(fwhm_xarr)
        if len(fwhm_yarr) > 2:
            self.fwhm_y = np.median(fwhm_yarr)
            if True:

                self.stats["stars"]["fwhm"] = (fwhm_xarr + fwhm_yarr)/2
                self.stats["stars"]["fwhm_xax"] = fwhm_xarr
                self.stats["stars"]["fwhm_yax"] = fwhm_yarr

                self.stats["fwhm"] = (self.fwhm_x + self.fwhm_y)/2.
                self.stats["fwhm_xax"] = self.fwhm_x
                self.stats["fwhm_yax"] = self.fwhm_y

        return self.fwhm_x, self.fwhm_y


    def sky_gradient(self,n_segments=10):
        image = self.image
        segments = self.make_segments(image)

        x = []
        y = []
        back = []
        back_mad = []
        for s in segments:
            x_tmp = (s["x"][0] + s["x"][1]) // 2
            y_tmp = (s["y"][0] + s["y"][1]) // 2
            x.append(x_tmp)
            y.append(y_tmp)
            back.append(np.median(s["subframe"]))
            back_mad.append(mad_std(s["subframe"]))

        x = np.array(x)
        y = np.array(y)
        back = np.array(back)

        A = np.vstack([np.ones_like(x), x, y, x ** 2, x * y, y ** 2]).T
        coeff, *_ = np.linalg.lstsq(A, back, rcond=None)

        x0 = min(x)
        y0 = min(y)
        xk = max(x)
        yk = max(y)

        le = []
        re = []
        ue = []
        de = []

        bk = []

        for xi, yi in zip(x, y):
            bk.append(self.polysurf(xi, yi, coeff))
            le.append(self.polysurf(x0, yi, coeff))
            re.append(self.polysurf(xk, yi, coeff))
            ue.append(self.polysurf(xi, yk, coeff))
            de.append(self.polysurf(xi, y0, coeff))

        max_amplitude = max(bk) - min(bk)
        max(le) - min(re)
        frame_gradient = max(max((max(le) - min(re)), (max(re) - min(le))),
                             max((max(ue) - min(de)), (max(de) - min(ue))))

        self.max_amplitude = max_amplitude
        self.frame_gradient = frame_gradient
        self.sky_surface_coeff = coeff
        self.sky_surface_bkg = bk
        self.sky_surface_x = x
        self.sky_surface_y = y

        self.stats["bkg_max_amplitude"] = self.max_amplitude
        self.stats["bkg_frame_gradient"] = self.frame_gradient
        self.stats["sky_surface_coeff"] = self.sky_surface_coeff

        tmp = {}
        tmp["surface_bkg"] = self.sky_surface_bkg
        tmp["surface_x"] = self.sky_surface_x
        tmp["surface_y"] = self.sky_surface_y
        self.stats["sky"] = tmp

    def hough_transform(self, maska, th_signal=100, steps=180):
        ys, xs = np.nonzero(maska)
        xs = xs.astype(float)
        ys = ys.astype(float)
        N = len(xs)

        ny, nx = self.image.shape
        rh0 = int((ny ** 2 + nx ** 2) ** 0.5)

        theta = np.deg2rad(np.linspace(-90, 90, steps))
        T = len(theta)

        cos_t = np.cos(theta)
        sin_t = np.sin(theta)

        rho_mtx = xs[:, None] * cos_t[None, :] + ys[:, None] * sin_t[None, :]

        rh = np.linspace(-rh0, rh0, 2 * rh0)
        ri = np.round(rho_mtx + rh0).astype(int)

        self.accumulator = np.zeros((len(rh), T))
        ri_flat = ri.ravel()
        ti_flat = np.tile(np.arange(T), N)
        np.add.at(self.accumulator, (ri_flat, ti_flat), 1)

        self.accumulator
        self.theta = theta
        self.rh = rh

        accu_max = self.accumulator.max(axis=1)
        peaks, _ = find_peaks(accu_max, height=th_signal)

        lines_rho = []
        lines_theta = []
        lines_val = []
        for p in peaks:
            rho = rh[p]
            ti = np.argmax(self.accumulator[p])
            t = self.theta[ti]
            val = self.accumulator[p, ti]
            lines_rho.append(rho)
            lines_theta.append(t)
            lines_val.append(val)

        idx = np.argsort(lines_val)[::-1]  # descending
        lines_rho = np.array(lines_rho)[idx]
        lines_theta = np.array(lines_theta)[idx]
        lines_val = np.array(lines_val)[idx]

        tmp = {}
        tmp["val"] = lines_val
        tmp["rho"] = lines_rho
        tmp["theta"] = lines_theta
        self.stats["lines"] = tmp


        return lines_val, lines_theta, lines_rho

    def line_filter(self,image,kernel1_size=3,kernel2_size=7,th1=3,th2=3):
        image = image

        kernel_l, kernel_r = self.line_detection_kernel(kernel1_size)  # tu sie zmienia
        result_l = convolve(image, kernel_l)
        maska_1l = result_l > np.median(result_l) + th1 * mad_std(result_l)
        result_r = convolve(image, kernel_r)
        maska_1r = result_r > np.median(result_r) + th1 * mad_std(result_r)

        kernel_l, kernel_r = self.line_detection_kernel(kernel2_size)  # tu sie zmienia
        result_l = convolve(image, kernel_l)
        maska_2l = result_l > np.median(result_l) + th2 * mad_std(result_r)
        result_r = convolve(image, kernel_r)
        maska_2r = result_r > np.median(result_r) + th2 * mad_std(result_r)

        maska = (maska_1l & maska_2l) ^ (maska_1r & maska_2r)

        return maska

    def make_segments(self,image, n_segments=10, overlap=0):
        height, width = image.shape
        seg_h = height // n_segments
        seg_w = width // n_segments

        segments = []

        for i in range(n_segments):
            for j in range(n_segments):
                result = {}

                y_start = (i * seg_h) - overlap
                if y_start < 0: y_start = 0
                x_start = (j * seg_w) - overlap
                if x_start < 0: x_start = 0
                y_end = ((i + 1) * seg_h) + overlap if i < n_segments - 1 else height
                x_end = ((j + 1) * seg_w) + overlap if j < n_segments - 1 else width

                subframe = image[y_start:y_end, x_start:x_end]

                result = {}
                result["x"] = [x_start, x_end]
                result["y"] = [y_start, y_end]
                result["subframe"] = subframe
                segments.append(result)

        return segments

    def polysurf(self,x, y, coeff):
        a0, a1, a2, a3, a4, a5 = coeff
        return a0 + a1 * x + a2 * y + a3 * x ** 2 + a4 * x * y + a5 * y ** 2

    def line_detection_kernel(self,kernel_half_size=5):
        size = 2 * kernel_half_size + 1
        center = kernel_half_size
        y, x = np.ogrid[:size, :size]
        distance = np.sqrt((x - center) ** 2 + (y - center) ** 2)
        inner = kernel_half_size - 1 / 2
        outer = kernel_half_size + 1 / 2

        left = ((x <= center) & (y <= center)) | ((x >= center) & (y >= center))
        right = ((x < center) & (y > center)) | ((x > center) & (y < center))

        # left
        mk = ((distance >= inner) & (distance <= outer))
        kernel_l = mk.astype(float)
        tmp_mk = ((distance >= inner) & (distance <= outer) & right)
        kernel_l[tmp_mk] = -1

        # right
        mk = ((distance >= inner) & (distance <= outer))
        kernel_r = mk.astype(float)
        tmp_mk = ((distance >= inner) & (distance <= outer) & left)
        kernel_r[tmp_mk] = -1

        return kernel_l, kernel_r

    def gauss_kernel(self, size, sigma):
        kernel = np.fromfunction(lambda x, y: (1 / (2 * np.pi * sigma ** 2)) * np.exp(
            -((x - (size - 1) / 2) ** 2 + (y - (size - 1) / 2) ** 2) / (2 * sigma ** 2)), (size, size))
        return kernel / np.sum(kernel)