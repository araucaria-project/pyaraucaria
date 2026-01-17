import numpy as np
from scipy.signal import find_peaks, fftconvolve
from scipy.ndimage.filters import maximum_filter
from scipy.ndimage import convolve
from scipy.ndimage import gaussian_filter

from astropy.stats import mad_std

class FFS:

    def __init__(self, image, gain=1., rn_noise=0.):
        self.image = np.transpose(image)
        self.gain = float(gain)
        self.rn_noise = float(rn_noise)
        self.saturation = 50000
        self.stats = {}


    def mk_stats(self):
        img = self.image.ravel()

        self.min = img.min()
        self.max = img.max()
        self.mean = img.mean()
        self.rms = img.std()

        s = np.sort(img)
        n = s.size
        q16 = s[int(0.159 * n)]
        q50 = s[n // 2]
        q84 = s[int(0.841 * n)]

        self.median = q50
        self.q_sigma_lower = q50 - q16
        self.q_sigma_upper = q84 - q50
        self.q_sigma = (q84 - q16) / 2.0

        self.sigma_quantile = self.q_sigma

        self.noise = np.sqrt(self.median / self.gain + self.rn_noise)

        self.stats = {
            "min": self.min,
            "max": self.max,
            "mean": self.mean,
            "median": self.median,
            "rms": self.rms,
            "q_sigma_lower": self.q_sigma_lower,
            "q_sigma_upper": self.q_sigma_upper,
            "q_sigma": self.q_sigma,
            "noise": self.noise,
        }

        self.stats_description = {
            "min": "Minimal pixel value in the image",
            "max": "Maximum pixel value in the image",
            "mean": "Mean pixel value (arithmetic average)",
            "median": "Median pixel value (robust background estimator)",
            "rms": "Standard deviation of pixel values",
            "q_sigma_lower": "Lower 1-sigma estimate from 15.9% quantile",
            "q_sigma_upper": "Upper 1-sigma estimate from 84.1% quantile",
            "q_sigma": "Robust sigma estimated from 15.9–84.1% quantiles",
            "noise": "Expected total noise (Poisson + read noise)",
        }

    def find_stars(self, threshold=5, method="sigma quantile", kernel_size=30, fwhm=10):

        self.coo = []
        self.adu = []
        self.threshold = float(threshold)
        self.method = method
        self.fwhm_adopted = float(fwhm)
        self.kernel_size = int(kernel_size)
        self.kernel_sigma = float(fwhm) / 2.355

        if self.method == "rms Poisson":
            self.sigma = self.noise
        elif self.method == "rms":
            self.sigma = self.rms
        elif self.method == "sigma quantile":
            self.sigma = self.q_sigma
        else:
            raise ValueError(f"Invalid method type {self.method}")

        mask1 = self.image > self.median + self.threshold * self.sigma
        data2 = gaussian_filter(self.image, sigma=self.kernel_sigma)
        mask2 = data2 == maximum_filter(data2, size=3)
        mask = mask1 & mask2

        coo = np.column_stack(np.nonzero(mask))
        val = self.image[mask]

        self.stats["stars"] = {}
        if len(coo) > 0:
            sorted_i = np.argsort(val)[::-1]
            self.coo = coo[sorted_i]
            self.adu = val[sorted_i]

            self.stats["stars"] = {
                "x": self.coo[:, 0],
                "y": self.coo[:, 1],
                "max_adu": self.adu
            }
        else:
            self.stats["stars"] = {
                "x": [],
                "y": [],
                "max_adu": []}

        self.stats_description["stars"] = {
        "x": "X coordinates of detected stars (pixel indices, 0-based)",
        "y": "Y coordinates of detected stars (pixel indices, 0-based)",
        "max_adu": "Peak ADU (brightness) of each detected star"
        }
        return self.coo, self.adu



    def fwhm(self, radius=10, all_stars=True):
        # deprecated
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
            if self.adu[i] < int(self.saturation) and i < i_max:
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

    def star_info(self,box=10,N_stars=None):

        self.ellipticity = []
        self.theta = []
        self.fwhm = []
        self.fwhm_x = []
        self.fwhm_y = []
        self.cpe = []

        if N_stars == None:
            N_stars = len(self.coo)

        ni = 0
        for i, (x, y) in enumerate(self.coo):
            if ni <= N_stars:

                e = np.nan
                t = np.nan
                f = np.nan
                fx = np.nan
                fy = np.nan
                cpe = np.nan

                if self.adu[i] < self.saturation:
                    ni = ni + 1

                    cut = self.image[x - box:x + box, y - box:y + box]

                    if cut.shape[0] > box-1 and cut.shape[1] > box-1:

                        cut = cut - np.median(cut[0, :])
                        _, e, t = ffs_pca(cut)

                        cut = cut - 0.5 * np.max(cut)

                        fx, fy = ffs_fwhm(cut)

                        if fx is not np.nan and fy is not np.nan:
                            f = (fx+fy)/2.
                            r = int(f)
                            cut = self.image[x - r:x + r, y - r:y + r]
                            cpe = ffs_cpe(cut)

                self.ellipticity.append(e)
                self.theta.append(t)
                self.fwhm.append(f)
                self.fwhm_x.append(fx)
                self.fwhm_y.append(fy)
                self.cpe.append(cpe)

            self.stats["stars"]["fwhm"] = self.fwhm
            self.stats["stars"]["fwhm_x"] = self.fwhm_x
            self.stats["stars"]["fwhm_y"] = self.fwhm_y
            self.stats["stars"]["ellipticity"] = self.ellipticity
            self.stats["stars"]["theta"] = self.theta
            self.stats["stars"]["cpe"] = self.cpe

            self.stats_description["stars"].update({

                "fwhm": (
                    "Full Width at Half Maximum (FWHM) of each detected star "
                ),

                "fwhm_x": (
                    "Full Width at Half Maximum (FWHM) of each detected star measured "
                    "along the X axis of the fitted profile (in pixels)"
                ),

                "fwhm_y": (
                    "Full Width at Half Maximum (FWHM) of each detected star measured "
                    "along the Y axis of the fitted profile (in pixels)"
                ),

                "ellipticity": (
                    "Ellipticity of each detected star"
                    "values close to 0 indicate round stars, "
                    "higher values indicate elongated profiles"
                ),

                "theta": (
                    "Position angle of the semi-major axis of each detected star, "
                    "measured counter-clockwise from the X axis of the detector "
                    "(in radians)"
                ),

                "cpe": (
                    "Central Pixel Excess (CPE) of each detected star, defined as "
                    "the contrast of the peak pixel relative to the local background "
                    "and background noise; higher values indicate sharper, more "
                    "centrally concentrated profiles"
                ),
            })

        return self.fwhm, self.fwhm_x, self.fwhm_y, self.ellipticity, self.theta, self.cpe


    def sky_gradient(self,n_segments=10):
        image = self.image
        segments = ffs_make_segments(image)

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
            bk.append(ffs_polysurf(xi, yi, coeff))
            le.append(ffs_polysurf(x0, yi, coeff))
            re.append(ffs_polysurf(xk, yi, coeff))
            ue.append(ffs_polysurf(xi, yk, coeff))
            de.append(ffs_polysurf(xi, y0, coeff))

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

        self.stats_description.update({

            "bkg_max_amplitude": (
                "Maximum amplitude of the background signal across the image, "
                "representing the highest estimated background level (in ADU)"
            ),

            "bkg_frame_gradient": (
                "Global background gradient across the image frame, describing "
                "large-scale background variation (e.g. due to sky brightness "
                "gradient, moonlight, or scattered light)"
            ),

            "sky_surface_coeff": (
                "Coefficients of the fitted background surface model describing "
                "the large-scale sky background variation across the image"
            ),

            "sky": {
                "surface_bkg": (
                    "Estimated sky background surface evaluated over the image grid, "
                    "representing the modeled large-scale background (in ADU)"
                ),
                "surface_x": (
                    "X-coordinate grid used for evaluating the sky background surface "
                    "(pixel indices, 0-based)"
                ),
                "surface_y": (
                    "Y-coordinate grid used for evaluating the sky background surface "
                    "(pixel indices, 0-based)"
                ),
            },
        })


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

        self.stats_description["lines"] = {

            "val": (
                "Detection strength or significance of each detected line, "
                "representing the response value of the line detection algorithm "
                "(higher values indicate more prominent linear features)"
            ),

            "rho": (
                "Perpendicular distance (ρ) of each detected line from the origin "
                "of the image coordinate system, as defined in the Hough transform "
                "parameter space (in pixels)"
            ),

            "theta": (
                "Orientation angle (θ) of each detected line, measured relative to "
                "the X axis of the image coordinate system (in radians)"
            ),
        }

        return lines_val, lines_theta, lines_rho


def ffs_pca(image_cut):
    e = np.nan
    t = np.nan
    f = np.nan
    lx, ly = image_cut.shape
    nx, ny = np.mgrid[0:lx, 0:ly]
    image_cut = image_cut.clip(min=0)
    I = image_cut.clip(min=0)
    It = I.sum()
    if It > 0:
        dx = nx - int(lx/2)
        dy = ny - int(ly/2)

        Mxx = np.sum(I * dx * dx) / It
        Myy = np.sum(I * dy * dy) / It
        Mxy = np.sum(I * dx * dy) / It

        cov = np.array([[Mxx, Mxy], [Mxy, Myy]])
        eigvals, eigvecs = np.linalg.eigh(cov)
        b2, a2 = eigvals

        if a2 > 0 and b2 > 0:

            f1 = 2.355 * np.sqrt(a2)
            f2 = 2.355 * np.sqrt(b2)
            f = (f1 + f2)/2

            e = 1.0 - np.sqrt(b2 / a2)
            vx, vy = eigvecs[:, 1]
            t = np.arctan2(vy, vx)
    return f,e,t

def ffs_cpe(image_cut):
    cpe = np.nan
    # central pixel excess
    cx = int(image_cut.shape[0]/2)
    cy = int(image_cut.shape[1]/2)

    if cx > 1 and cy > 1:
        x_max, y_max = np.unravel_index(np.argmax(image_cut), image_cut.shape)
        I_max = image_cut[x_max, y_max]

        bck = image_cut.copy().astype(float)
        bck[x_max,y_max] = np.nan

        I_bck = np.nanmean(image_cut)
        I_std = np.nanstd(image_cut)

        # Zabezpieczenie przed dzieleniem przez zero
        if I_std == 0 or np.isnan(I_std):
            return np.nan

        cpe = (I_max - I_bck) / I_std

    return cpe

def ffs_fwhm(image_cut):
    cx = int(image_cut.shape[0]/2)
    cy = int(image_cut.shape[1]/2)
    line_x = image_cut[cx, :]
    line_y = image_cut[:, cy]

    fwhm_x = ffs_fwhm_1d(line_x)
    fwhm_y = ffs_fwhm_1d(line_y)

    return fwhm_x,fwhm_y



def ffs_fwhm_1d(line):

    idx = np.where(line > 0)[0]
    if len(idx) < 2:
        return np.nan

    i1 = idx[0]
    i2 = idx[-1]

    # interpolacja lewego brzegu
    if i1 > 0:
        x1, x2 = i1-1, i1
        y1, y2 = line[x1], line[x2]
        xl = x1 + (0 - y1) / (y2 - y1)
    else:
        xl = i1

    # interpolacja prawego brzegu
    if i2 < len(line)-1:
        x1, x2 = i2, i2+1
        y1, y2 = line[x1], line[x2]
        xr = x1 + (0 - y1) / (y2 - y1)
    else:
        xr = i2

    return xr - xl



def ffs_make_segments(image, n_segments=10, overlap=0):
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

def ffs_polysurf(x, y, coeff):
    a0, a1, a2, a3, a4, a5 = coeff
    return a0 + a1 * x + a2 * y + a3 * x ** 2 + a4 * x * y + a5 * y ** 2

def ffs_line_filter(image,kernel1_size=3,kernel2_size=7,th1=3,th2=3):
    image = image

    kernel_l, kernel_r = ffs_line_detection_kernel(kernel1_size)  # tu sie zmienia
    result_l = convolve(image, kernel_l)
    maska_1l = result_l > np.median(result_l) + th1 * mad_std(result_l)
    result_r = convolve(image, kernel_r)
    maska_1r = result_r > np.median(result_r) + th1 * mad_std(result_r)

    kernel_l, kernel_r = ffs_line_detection_kernel(kernel2_size)  # tu sie zmienia
    result_l = convolve(image, kernel_l)
    maska_2l = result_l > np.median(result_l) + th2 * mad_std(result_r)
    result_r = convolve(image, kernel_r)
    maska_2r = result_r > np.median(result_r) + th2 * mad_std(result_r)

    maska = (maska_1l & maska_2l) ^ (maska_1r & maska_2r)

    return maska

def ffs_line_detection_kernel(kernel_half_size=5):
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

def ffs_gauss2d_kernel(size, sigma):
    kernel = np.fromfunction(lambda x, y: (1 / (2 * np.pi * sigma ** 2)) * np.exp(
        -((x - (size - 1) / 2) ** 2 + (y - (size - 1) / 2) ** 2) / (2 * sigma ** 2)), (size, size))
    return kernel / np.sum(kernel)