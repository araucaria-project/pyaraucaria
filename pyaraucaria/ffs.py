import warnings

import numpy as np
from scipy.ndimage import maximum_filter, median_filter
from scipy.ndimage import convolve
from scipy.ndimage import gaussian_filter
from scipy.ndimage import binary_closing, label, map_coordinates
from scipy.interpolate import LinearNDInterpolator, NearestNDInterpolator
from scipy.special import erf

from astropy.stats import mad_std
from astropy.io import fits
from astropy.table import Table



class FFS:
    """Fast FITS Statistics: quick frame statistics, sky background, star and trail detection, masks.

    FFS gives estimates (quality control, masks for photometry); it does not replace photutils.

    Parameters
    ----------
    image : 2D array
        Image in ADU.
    gain : float
        Gain in e-/ADU. Used by find_stars(method="rms Poisson").
    rn_noise : float
        Read-out noise. Used by find_stars(method="rms Poisson").
    saturation : float
        Saturation level in ADU. Pixels at or above it are saturated (mk_saturation_mask);
        diffraction spikes are recognised only on saturated stars (classify_lines).
    bin : int
        Internal bin x bin binning (mean). 1 = no binning. Ignored when autobin is True.
    autobin : bool
        Choose bin automatically so that the longer side of the binned image is
        <= FFS.AUTOBIN_LIMIT (2500 px), e.g. 4096 x 4096 -> bin 2. Recommended for large frames.

    Notes
    -----
    Binning and memory: by default (bin=1, autobin=False) the image is NOT binned. On large
    frames (e.g. 4096 x 4096) find_lines() - also called by calc_frame_fwhm() - then needs
    more than 10 GB of RAM. Use autobin=True (or bin >= 2) for frames larger than ~2500 px.
    Binning is the user's decision - FFS never switches it on by itself.

    Computed on the binned image (self.binned_image): sky map (sky_map), threshold mask,
    lines (find_lines, hough_transform) and star detection (find_stars). Computed on the
    original: saturation, cosmics (mk_cosmics_mask), star measurements (star_info: FWHM,
    ellipticity), classify_lines, calibration maps. Everything returned to the user (star
    and line coordinates, FWHM, masks in self.masks, maps in self.maps) is in the frame and
    size of the original, and pixel parameters (fwhm, min_length, max_gap, box, ...) are
    given in original pixels. Exceptions: find_lines(half_width, max_shift) are in binned
    pixels. Internal, binned size: self.binned_image, self.binned_masks,
    self.accumulator / self.rh / self.hough_theta.

    Masks and maps: self.masks holds boolean masks of the image size, one per method
    (described in self.masks_description). Each method takes exclude - the names of the masks
    whose pixels it skips (bad_mask()); masks that do not exist yet are ignored. The defaults
    are built from FFS.BAD_PIXELS and FFS.SOURCES, e.g. a stars mask is skipped in statistics
    and background but not in star detection. self.maps holds maps (e.g. "sky") of the image
    size, described in self.maps_description.

    Background sigma in find_stars(method=...): "sigma quantile" (default) - mad_std of the
    residual image - sky (_background); "rms" - standard deviation of good pixels (mk_stats);
    "rms Poisson" - sqrt(median / gain + rn_noise^2) (mk_stats).
    """

    AUTOBIN_LIMIT = 2500   # autobin: najmniejszy bin, przy ktorym dluzszy bok obrazu zbinowanego <= limit

    # save_mef: krotkie (<= 8 znakow) klucze FITS dla skalarow z stats["frame"]
    FITS_KEYS = {
        "min": "FFMIN", "max": "FFMAX", "mean": "FFMEAN", "median": "FFMEDIAN", "rms": "FFRMS",
        "q_sigma_lower": "FFQSIGL", "q_sigma_upper": "FFQSIGU", "q_sigma": "FFQSIG",
        "noise": "FFNOISE", "n_bad": "FFNBAD",
        "fwhm": "FFFWHM", "fwhm_x": "FFFWHMX", "fwhm_y": "FFFWHMY", "ellipticity": "FFELLIP",
        "theta_spread": "FFTHSPR", "cpe": "FFCPE", "shape": "FFSHAPE", "ci": "FFCI",
        "used_stars": "FFNSTARS",
        "bkg_max_amplitude": "FFBKGAMP", "bkg_frame_gradient": "FFBKGGRD",
        "psf_e1": "FFPSFE1", "psf_e2": "FFPSFE2", "psf_e_mean": "FFPSFE", "psf_theta_mean": "FFPSFTH",
        "psf_e_t": "FFPSFET", "psf_e_x": "FFPSFEX",
        "psf_fwhm_tilt": "FFFWTILT", "psf_fwhm_tilt_angle": "FFFWTANG", "psf_fwhm_curv": "FFFWCURV",
    }
    # klucze strukturalne naglowka oryginalu, nie kopiowane do PRIMARY (NAXISn osobno)
    FITS_STRUCTURAL = {"SIMPLE", "BITPIX", "EXTEND", "BZERO", "BSCALE", "PCOUNT", "GCOUNT",
                       "XTENSION", "EXTNAME", "CHECKSUM", "DATASUM"}

    # nazwy masek do domyslnych list exclude: zle piksele i zrodla (slady, spajki, gwiazdy)
    BAD_PIXELS = ("nonfinite", "saturation", "cosmics", "hot_pixels", "cold_pixels",
                  "bad_columns", "bad_rows")
    SOURCES = ("lines", "spikes", "stars")
    EXCLUDE_ALL = BAD_PIXELS + SOURCES
    # detekcja gwiazd: bez saturacji i gwiazd (nasycone gwiazdy maja zostac w coo)
    EXCLUDE_DETECTION = ("nonfinite", "cosmics", "hot_pixels", "cold_pixels",
                         "bad_columns", "bad_rows", "lines", "spikes")
    # pomiar gwiazd: bez gwiazd
    EXCLUDE_MEASURE = BAD_PIXELS + ("lines", "spikes")
    # klasyfikacja linii: bez saturacji i linii (to je badamy) i gwiazd
    EXCLUDE_LINES = ("nonfinite", "cosmics", "hot_pixels", "cold_pixels", "bad_columns", "bad_rows")

    def __init__(self, image, gain=1., rn_noise=0., saturation = 50000, bin=1, autobin=False):
        self.image = image
        self.gain = float(gain)
        self.rn_noise = float(rn_noise)
        self.saturation = saturation

        # binning tylko do uzytku wewnetrznego (tlo, detekcja, linie); wyniki dla uzytkownika
        # sa w wymiarach oryginalu. autobin=True dobiera bin sam, ignorujac argument bin.
        if autobin:
            bin = int(np.ceil(max(np.shape(image)) / FFS.AUTOBIN_LIMIT))
        self.bin = max(1, int(bin))
        self.binned_image = FFS.bin_image(image, self.bin)

        self.min = None
        self.max = None
        self.mean = None
        self.rms = None
        self.median = None
        self.q_sigma = None
        self.q_sigma_lower = None
        self.q_sigma_upper = None
        self.sigma_quantile = None
        self.noise = None

        self.stars = None
        self.coo = None
        self.adu = None

        self.fs_threshold = None
        self.fs_method = None
        self.fs_fwhm_adopted = None
        self.fs_kernel_sigma = None

        self.box_mag = None
        self.bkg = None
        self.ellipticity = None
        self.theta = None
        self.fwhm = None
        self.fwhm_x = None
        self.fwhm_y = None
        self.cpe = None

        self.max_amplitude = None
        self.frame_gradient = None
        self.sky_surface_coeff = None

        self.sky_surface_bkg = None
        self.sky_surface_x = None
        self.sky_surface_y = None
        self.sky = None

        self.lines_val = None
        self.lines_theta = None
        self.lines_rho = None

        self.frame_fwhm = None
        self.frame_fwhm_x = None
        self.frame_fwhm_y = None
        self.frame_ellipticity = None
        self.frame_theta_spread = None
        self.frame_cpe = None
        self.frame_shape = None
        self.frame_ci = None

        self.rh = None
        self.stats = {
            "frame": {},
            "stars": {}
        }
        self.stats_description = {
            "frame": {},
            "stars": {}
        }

        self.lines = None              # tabela wykrytych linii (find_lines)
        self.psf_grid = None           # tabela PSF w wezlach rzadkiej siatki - "wasy" (mk_star_maps)
        self.masks = {}                # maski logiczne w rozmiarze obrazu, osobno dla kazdej metody
        self.binned_masks = {}         # maski wewnetrzne w rozmiarze binned_image (np. wejscie do Hougha)
        self.masks_description = {}

        nonfinite = ~np.isfinite(self.image)
        if nonfinite.any():            # NaN/inf na wejsciu (np. po redukcji) traktujemy jak zle piksele
            self.masks["nonfinite"] = nonfinite
            self.masks_description["nonfinite"] = "Non-finite input pixels (NaN or inf)"

        self.maps = {}
        self.maps_description = {}

    def bad_mask(self, exclude=EXCLUDE_ALL):
        """Union (OR) of the masks named in exclude; names of masks that do not exist are ignored."""
        bad = np.zeros(self.image.shape, dtype=bool)
        for name in exclude:
            if name in self.masks:
                bad |= self.masks[name]
        return bad

    def save_mef(self, path, header=None, overwrite=False):
        """Save maps, masks and tables to one multi-extension FITS file.

        PRIMARY: no data; a copy of header (if given, e.g. the original frame header),
        FFS parameters and the scalar frame statistics under short FF* keys (FFS.FITS_KEYS).
        Then one float32 ImageHDU per map (EXTNAME = map name in upper case), MASKS - int32
        bit mask (bit i = mask MASKi from its header) and the tables STARS, LINES, SKYSURF (sky_gradient),
        PSFGRID (mk_star_maps). Missing products are skipped.
        """
        primary = fits.PrimaryHDU()
        if header is not None:
            # bez kluczy strukturalnych oryginalu (NAXISn, BZERO...) - PRIMARY nie ma danych
            for card in header.cards:
                if card.keyword not in FFS.FITS_STRUCTURAL and not card.keyword.startswith("NAXIS"):
                    primary.header.append(card)

        params = {"FFBIN": (self.bin, "FFS internal binning"),
                  "FFGAIN": (self.gain, "FFS gain [e-/ADU]"),
                  "FFRN": (self.rn_noise, "FFS read-out noise"),
                  "FFSATUR": (self.saturation, "FFS saturation level [ADU]")}
        for key, (value, comment) in params.items():
            primary.header[key] = (value, comment)

        # tylko skonczone skalary - NaN/inf nie sa dozwolone w naglowku FITS
        descr = self.stats_description["frame"]
        for name, value in self.stats["frame"].items():
            if name in FFS.FITS_KEYS and np.isscalar(value) and np.isfinite(value):
                primary.header[FFS.FITS_KEYS[name]] = (value, FFS._fits_text(descr.get(name, ""))[:46])

        hdus = [primary]
        for name, data in self.maps.items():
            hdu = fits.ImageHDU(np.asarray(data, dtype=np.float32), name=name.upper())
            hdu.header["COMMENT"] = FFS._fits_text(self.maps_description.get(name, ""))
            hdus.append(hdu)

        if self.masks:
            if len(self.masks) > 31:
                raise ValueError(f"save_mef: {len(self.masks)} masks do not fit in int32 (max 31)")
            bits = np.zeros(self.image.shape, dtype=np.int32)
            hdu = fits.ImageHDU(name="MASKS")
            for i, (name, mask) in enumerate(self.masks.items()):
                bits |= mask.astype(np.int32) << i
                hdu.header[f"MASK{i}"] = (name, FFS._fits_text(self.masks_description.get(name, ""))[:46])
            hdu.data = bits
            hdus.append(hdu)

        for name, table in (("STARS", self.stars), ("LINES", self.lines), ("SKYSURF", self.sky),
                            ("PSFGRID", self.psf_grid)):
            if table is not None:
                hdu = fits.table_to_hdu(Table(table))
                hdu.name = name
                hdus.append(hdu)

        fits.HDUList(hdus).writeto(path, overwrite=overwrite)

    def _background(self, binned=False, exclude=EXCLUDE_ALL):
        """Poziom tla i sigma: (maps["sky"], sigma residuum). Mapa tla liczona, jesli jej jeszcze nie ma.

        Sigma to mad_std residuum image - sky dobrych pikseli: q_sigma calej klatki zawiera tez
        rozrzut samego gradientu tla, a MAD jest odporniejszy na gwiazdy niz kwantyle 15.9/84.1
        (gesty kadr: kwantyle zawyzaja sigma kilkukrotnie).
        binned=True: tlo i sigma dla binned_image (zbinowana mapa tla; sigma jest mniejsza,
        bo binning usrednia szum).
        """
        if "sky" not in self.maps:
            self.sky_map()
        if binned:
            image, sky = self.binned_image, FFS.bin_image(self.maps["sky"], self.bin)
            bad = self.to_binned_mask(self.bad_mask(exclude))
        else:
            image, sky, bad = self.image, self.maps["sky"], self.bad_mask(exclude)
        residual = image - sky
        good = residual[~bad] if bad.sum() < bad.size else residual.ravel()
        return sky, mad_std(good), "sky map"

    def mk_stats(self, exclude=EXCLUDE_ALL):
        """Basic frame statistics -> self.stats["frame"] (min, max, mean, median, rms, quantile sigmas, noise).

        Computed from good pixels (outside bad_mask(exclude)); min and max from all finite pixels,
        so that saturation stays visible.

        Parameters
        ----------
        exclude : sequence of str
            Names of the masks whose pixels are skipped (bad_mask); missing masks are ignored.
        """
        # statystyki z dobrych pikseli (poza bad_mask); min/max ze wszystkich, zeby bylo widac saturacje
        bad = self.bad_mask(exclude)
        self.n_bad = int(bad.sum())
        img = self.image[~bad] if self.n_bad < bad.size else self.image.ravel()

        finite = self.image[~self.masks["nonfinite"]] if "nonfinite" in self.masks else self.image
        self.min = finite.min()
        self.max = finite.max()
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

        self.noise = np.sqrt(self.median / self.gain + self.rn_noise**2)

        self.stats["frame"] = {
            "min": self.min,
            "max": self.max,
            "mean": self.mean,
            "median": self.median,
            "rms": self.rms,
            "q_sigma_lower": self.q_sigma_lower,
            "q_sigma_upper": self.q_sigma_upper,
            "q_sigma": self.q_sigma,
            "noise": self.noise,
            "n_bad": self.n_bad,
        }

        self.stats_description["frame"] = {
            "min": "Minimal pixel value in the image",
            "max": "Maximum pixel value in the image",
            "mean": "Mean pixel value (arithmetic average)",
            "median": "Median pixel value (robust background estimator)",
            "rms": "Standard deviation of pixel values",
            "q_sigma_lower": "Lower 1-sigma estimate from 50% - 15.9% quantile",
            "q_sigma_upper": "Upper 1-sigma estimate from 84.1% - 50% quantile",
            "q_sigma": "Robust sigma estimated from 15.9–84.1% quantiles",
            "noise": "Expected total noise (Poisson + read noise)",
            "n_bad": "Number of excluded pixels (bad_mask); other statistics except min/max use only the remaining pixels",
        }

    def mk_threshold_mask(self, nsigma=3.0, exclude=EXCLUDE_ALL):
        """Mask of pixels above the background: image > sky + nsigma * sigma -> masks["threshold"].

        The sky is self.maps["sky"] (computed if missing), sigma is mad_std of the residual
        image - sky. Computed on the binned image (binned_masks["threshold"] is the input of
        the Hough transform in find_lines); masks["threshold"] is the same mask in original size.

        Parameters
        ----------
        nsigma : float
            Threshold above the sky map in units of background sigma. Lower values let
            fainter trails into line detection, but also more noise.
        exclude : sequence of str
            Names of the masks whose pixels are skipped (bad_mask); missing masks are ignored.
        """
        level, sigma, level_name = self._background(binned=True, exclude=exclude)
        bad = self.to_binned_mask(self.bad_mask(exclude))
        self.binned_masks["threshold"] = (self.binned_image > level + nsigma * sigma) & ~bad
        self.masks["threshold"] = self.to_original_mask(self.binned_masks["threshold"])
        self.masks_description["threshold"] = (
            f"Pixels above {level_name} + {nsigma} * sigma (sigma = {sigma:.2f} ADU, robust; "
            f"computed on the image binned {self.bin}x{self.bin}; "
            f"excluding bad_mask; input mask for line detection)"
        )

    def mk_saturation_mask(self, saturation=None):
        """Mask of saturated pixels: image >= saturation -> masks["saturation"] (excluding mask).

        Parameters
        ----------
        saturation : float | None
            Saturation level in ADU. None uses self.saturation (constructor).
        """
        if saturation is None:
            saturation = self.saturation
        self.masks["saturation"] = self.image >= saturation
        self.masks_description["saturation"] = f"Saturated pixels (value >= {saturation} ADU)"

    def mk_columns_map(self, block=64, window=15, exclude=EXCLUDE_ALL):
        """Maps of column and row deviations -> maps["columns"], maps["rows"] (calibration frames:
        bias, dark, flat, no stars).

        The sky map is subtracted first (sky_map, computed if missing; the number of segments
        can be set before with sky_map(n_segments=...)). maps["columns"] is the median of a
        column in a block of rows minus the median of the neighbouring columns, in ADU, image
        size; maps["rows"] the same for rows. Thresholded masks: mk_bad_columns_mask().

        Parameters
        ----------
        block : int
            Length (in pixels along the column) of the blocks in which column medians are computed.
        window : int
            Number of neighbouring columns whose median is the reference.
        exclude : sequence of str
            Names of the masks whose pixels are skipped (bad_mask); missing masks are ignored.
        """
        if "sky" not in self.maps:
            self.sky_map()
        resid = np.where(self.bad_mask(exclude), np.nan, self.image - self.maps["sky"])
        self.maps["columns"] = FFS._column_deviation(resid, block, window)
        self.maps["rows"] = FFS._column_deviation(resid.T, block, window).T
        for name, what in (("columns", "column"), ("rows", "row")):
            self.maps_description[name] = (
                f"{what.capitalize()} median in blocks of {block} px minus the median of {window} neighbouring "
                f"{what}s, sky map subtracted (ADU)"
            )

    def mk_bad_columns_mask(self, threshold):
        """Masks of bad columns and rows -> masks["bad_columns"], masks["bad_rows"] (excluding masks).

        Requires mk_columns_map().

        Parameters
        ----------
        threshold : float
            A pixel is bad when |maps["columns"]| (or |maps["rows"]|) exceeds this value, in ADU.
        """
        for name, what in (("columns", "column"), ("rows", "row")):
            self.masks["bad_" + name] = np.abs(self.maps[name]) > threshold
            self.masks_description["bad_" + name] = (
                f"Bad {what}s: |maps['{name}']| > {threshold} ADU (mk_bad_columns_mask)"
            )

    def mk_pixels_map(self, size=5):
        """Map of single-pixel deviations: image - median_filter(image, size x size) -> maps["pixels"], in ADU.

        For calibration frames (preferably a master; on a single dark a hot pixel and a cosmic
        ray look the same). Thresholded masks: mk_hot_pixels_mask(), mk_cold_pixels_mask().

        Parameters
        ----------
        size : int
            Size of the median filter box, in pixels.
        """
        image = self.image.astype(float)
        self.maps["pixels"] = image - median_filter(image, size=size)
        self.maps_description["pixels"] = f"Pixel value minus the median of its {size}x{size} neighbourhood (ADU)"

    def mk_hot_pixels_mask(self, threshold):
        """Mask of hot pixels: maps["pixels"] > threshold -> masks["hot_pixels"] (excluding mask).

        Requires mk_pixels_map().

        Parameters
        ----------
        threshold : float
            Minimal excess over the local median, in ADU.
        """
        self.masks["hot_pixels"] = self.maps["pixels"] > threshold
        self.masks_description["hot_pixels"] = f"Hot pixels: maps['pixels'] > {threshold} ADU (mk_hot_pixels_mask)"

    def mk_cold_pixels_mask(self, threshold):
        """Mask of cold pixels: maps["pixels"] < -threshold -> masks["cold_pixels"] (excluding mask).

        Requires mk_pixels_map().

        Parameters
        ----------
        threshold : float
            Minimal deficit below the local median, in ADU.
        """
        self.masks["cold_pixels"] = self.maps["pixels"] < -threshold
        self.masks_description["cold_pixels"] = f"Cold pixels: maps['pixels'] < -{threshold} ADU (mk_cold_pixels_mask)"

    def _refine_on_original(self, coo, exclude=EXCLUDE_DETECTION):
        """Pozycje (y, x) z obrazu zbinowanego -> najjasniejszy dobry piksel oryginalu w otoczeniu
        bloku (3x3 bloki). Zwraca (coo, adu) w oryginale; adu to wartosc piksela oryginalu
        (srednia bloku nasyconej gwiazdy spada ponizej saturacji). Duplikaty usuwane, kolejnosc
        zachowana."""
        b = self.bin
        ny, nx = self.image.shape
        bad = self.bad_mask(exclude)
        out = np.empty_like(coo)
        for i, (yb, xb) in enumerate(coo):
            y0, y1 = max(0, (yb - 1) * b), min(ny, (yb + 2) * b)
            x0, x1 = max(0, (xb - 1) * b), min(nx, (xb + 2) * b)
            cut = np.where(bad[y0:y1, x0:x1], -np.inf, self.image[y0:y1, x0:x1])
            dy, dx = np.unravel_index(np.nanargmax(cut), cut.shape)
            out[i] = (y0 + dy, x0 + dx)
        _, first = np.unique(out, axis=0, return_index=True)
        out = out[np.sort(first)]
        return out, self.image[out[:, 0], out[:, 1]]

    def mk_cosmics_mask(self, nsig=5.0, contrast=3.0):
        """Mask of cosmic rays (single pixels sharper than the PSF) -> masks["cosmics"] (excluding mask).

        For science frames. A pixel is a cosmic ray when:
        1. image - median3x3 > nsig * sigma (it sticks out above its neighbours beyond noise),
        2. image - sky > contrast * (median3x3 - sky) (sharper than a star: for a star with
           FWHM >= 2.5 px the median of the neighbours is about half of the peak, ratio ~2).
        The denominator in 2. is bounded from below by sigma. Sky and sigma from the sky map.
        Stars with FWHM < ~2.5 px may be partly masked; cosmic-ray tracks wider than 1-2 px
        are caught only partly.

        Parameters
        ----------
        nsig : float
            Minimal excess over the 3x3 median, in units of background sigma.
        contrast : float
            Minimal ratio of the pixel to its 3x3 median (both above the sky); a star gives ~2.
        """
        level, sigma, _ = self._background()
        image = self.image.astype(float)
        med3 = median_filter(image, size=3)
        self.masks["cosmics"] = (
            ((image - med3) > nsig * sigma)
            & ((image - level) > contrast * np.maximum(med3 - level, sigma))
        )
        self.masks_description["cosmics"] = (
            f"Cosmic rays: image - median3x3 > {nsig} sigma and image - sky > {contrast} * (median3x3 - sky) "
            f"(mk_cosmics_mask)"
        )

    def mk_stars_mask(self, k=3.0, fwhm=None):
        """Mask of stars: discs of radius k * fwhm around all detections in self.coo -> masks["stars"].

        For statistics and background (in FFS.SOURCES); not used in star detection and
        measurement. Requires find_stars(); all detections are masked, also saturated and
        not measured ones (star_info keeps them with NaN).

        Parameters
        ----------
        k : float
            Disc radius in units of fwhm.
        fwhm : float | None
            FWHM in pixels; None = self.frame_fwhm (calc_star_stats / calc_frame_fwhm).
        """
        if fwhm is None:
            fwhm = self.frame_fwhm
        if fwhm is None or not np.isfinite(fwhm):
            raise ValueError("mk_stars_mask: no fwhm (run calc_frame_fwhm or give fwhm)")
        r = k * fwhm
        ri = int(np.ceil(r))
        ny, nx = self.image.shape
        mask = np.zeros((ny, nx), dtype=bool)
        for y, x in np.asarray(self.coo, dtype=int).reshape(-1, 2):
            y0, y1 = max(0, y - ri), min(ny, y + ri + 1)
            x0, x1 = max(0, x - ri), min(nx, x + ri + 1)
            yy, xx = np.ogrid[y0:y1, x0:x1]
            mask[y0:y1, x0:x1] |= (yy - y) ** 2 + (xx - x) ** 2 <= r ** 2
        self.masks["stars"] = mask
        self.masks_description["stars"] = (
            f"Stars: discs of radius {k} * FWHM = {r:.1f} px around {len(self.coo)} detections (mk_stars_mask)"
        )

    def find_stars(self, threshold=5, method="sigma quantile", fwhm=10,
                   min_smoothed_sigma=None, rank_by="raw",
                   max_concentration=None, aperture_radius=None,
                   min_pixels_above_threshold=None, exclude=EXCLUDE_DETECTION):
        """Detect star peaks. Detection runs on ``self.binned_image`` (= ``self.image``
        for bin=1); the returned ``self.coo`` (y, x) and ``self.adu`` are in the original
        image: each detection is moved to the brightest good original pixel near its block.

        Parameters
        ----------
        threshold : float
            Detection threshold in units of background sigma (raw image).
        method : {"sigma quantile", "rms", "rms Poisson"}
            How background sigma is estimated: "sigma quantile" - mad_std of image - sky map;
            "rms" - standard deviation of good pixels; "rms Poisson" - sqrt(median / gain + rn_noise^2).
        fwhm : float
            Adopted PSF FWHM (original pixels). Drives the Gaussian kernel sigma
            used to smooth the image before the local-maximum test.
        min_smoothed_sigma : float | None
            When set, also require the *smoothed-image* peak to exceed
            ``median(data2) + min_smoothed_sigma * mad(data2)``. This is
            a kernel-matched (matched-filter) SNR test that rejects hot
            pixels and narrow noise spikes — they collapse to a small
            fraction of their raw amplitude under the Gaussian
            convolution while real PSF-shaped sources stay bright.
            Default ``None`` preserves prior behaviour.
        max_concentration : float | None
            Reject candidates where the central pixel carries more than
            this fraction of the background-subtracted aperture flux.
            A hot pixel / cosmic ray / readout spike has concentration
            ≈ 1.0 (all signal in one pixel); a real PSF has
            concentration ≈ 0.2-0.4 depending on FWHM vs aperture.
            Default ``None`` disables. Recommended ~0.5 for typical
            guiders — cuts single-pixel artifacts regardless of how
            bright they are. Aperture from ``aperture_radius``.
        aperture_radius : int | None
            Half-width in original pixels of the square aperture used for both
            ``max_concentration`` and ``rank_by="aperture"``. ``None``
            picks ``max(2, round(2 * fwhm))`` — wide enough to capture
            most of an in-focus PSF.
        rank_by : {"raw", "smoothed", "aperture"}
            Sort key for the returned candidate list (descending). All
            three return the same coords + raw-peak ``adu`` array; only
            order differs. ``"raw"`` is the FFS legacy default — peak
            pixel ADU, hot pixels win. ``"smoothed"`` uses the
            kernel-matched amplitude (real stars outrank surviving hot
            pixels for bright targets, but dim real stars can still
            lose). ``"aperture"`` ranks by background-subtracted flux
            inside ``aperture_radius`` — most robust against
            single-pixel artifacts.
        min_pixels_above_threshold : int | None
            Reject candidates whose ``mask1`` (raw threshold) footprint
            inside the aperture box has fewer than this many pixels.
            Physically: a real PSF spans many pixels above the noise
            floor (e.g., FWHM=3 PSF at 5σ gives ~5-15 above-threshold
            pixels); a hot pixel, cosmic-ray streak, or random noise
            spike has 1-3. Cheapest possible "shape" filter — uses
            the already-computed mask1, just counts pixels per box.
            Default ``None`` disables. Recommended ~5 for typical
            in-focus PSFs.
        exclude : sequence of str
            Names of the masks whose pixels are skipped (bad_mask); missing masks are ignored.
            Without "saturation" and "stars" by default, so that saturated and already
            masked stars are detected.
        """

        self.coo = []
        self.adu = []
        self.fs_threshold = float(threshold) # inna nazwa zwiazana z instancja
        self.fs_method = method
        self.fs_fwhm_adopted = float(fwhm)
        self.fs_kernel_sigma = float(fwhm) / 2.355
        self.fs_min_smoothed_sigma = (
            float(min_smoothed_sigma) if min_smoothed_sigma is not None else None
        )
        self.fs_max_concentration = (
            float(max_concentration) if max_concentration is not None else None
        )
        self.fs_aperture_radius = (
            int(aperture_radius) if aperture_radius is not None
            else max(2, int(round(2.0 * float(fwhm))))
        )
        self.fs_rank_by = rank_by

        # detekcja na obrazie zbinowanym: fwhm i aperture_radius (podane w pikselach oryginalu,
        # tak tez zapisane w fs_*) przeliczone na piksele zbinowane; wynikowe coo/adu w oryginale
        b = self.bin
        kernel_sigma = self.fs_kernel_sigma / b
        ap_radius = max(1, int(round(self.fs_aperture_radius / b)))

        # tlo i sigma obrazu zbinowanego (mapa tla, mad_std residuum)
        level, q_sigma, _ = self._background(binned=True)

        # rms/noise z mk_stats() sa dla oryginalu; binning srednia zmniejsza szum ~ bin razy
        if self.fs_method == "rms Poisson":
            self.fs_sigma = self.noise / b
        elif self.fs_method == "rms":
            self.fs_sigma = self.rms / b
        elif self.fs_method == "sigma quantile":
            self.fs_sigma = q_sigma
        else:
            raise ValueError(f"Invalid method type {self.fs_method}")

        # Cast to float first: gaussian_filter preserves the input dtype, so on
        # integer FITS data (uint16/int16) the smoothed image — and every test
        # built on it (min_smoothed_sigma, rank_by="smoothed") — would be
        # truncated to integers.
        image = self.binned_image.astype(float)
        # zle piksele (bad_mask, NaN) wypelniamy tlem, zeby nie byly kandydatami i nie psuly
        # wygladzania; saturacji nie, bo nasycone gwiazdy maja zostac w coo/adu
        bad = self.to_binned_mask(self.bad_mask(exclude))
        image[bad] = level[bad] if np.ndim(level) else level

        mask1 = image > level + self.fs_threshold * self.fs_sigma
        data2 = gaussian_filter(image, sigma=kernel_sigma)
        mask2 = data2 == maximum_filter(data2, size=3)
        mask = mask1 & mask2

        if self.fs_min_smoothed_sigma is not None:
            # MAD-based robust sigma of the smoothed image; falls back to
            # std() if MAD collapses (uniform smoothed frame).
            smoothed_med = float(np.median(data2))
            smoothed_mad = float(np.median(np.abs(data2 - smoothed_med)))
            smoothed_sigma = 1.4826 * smoothed_mad
            if smoothed_sigma <= 0:
                smoothed_sigma = float(data2.std()) or 1.0
            mask3 = data2 > smoothed_med + self.fs_min_smoothed_sigma * smoothed_sigma
            mask = mask & mask3

        coo = np.column_stack(np.nonzero(mask))
        val = self.binned_image[mask]

        # Per-candidate aperture statistics (background-subtracted) —
        # needed for max_concentration filter and rank_by="aperture".
        # Computed once for both, since both are cheap (~50 ops/cand).
        H, W = image.shape
        n_cand = len(coo)
        aperture_excess = np.zeros(n_cand, dtype=float)
        peak_excess = np.zeros(n_cand, dtype=float)
        for i in range(n_cand):
            y, x = int(coo[i, 0]), int(coo[i, 1])
            y0 = max(0, y - ap_radius)
            y1 = min(H, y + ap_radius + 1)
            x0 = max(0, x - ap_radius)
            x1 = min(W, x + ap_radius + 1)
            bg_patch = level[y0:y1, x0:x1] if np.ndim(level) else level
            bg_peak = level[y, x] if np.ndim(level) else level
            aperture_excess[i] = float(np.sum(image[y0:y1, x0:x1] - bg_patch))
            peak_excess[i] = float(image[y, x] - bg_peak)

        # Concentration-index cull — single-pixel artifacts have all
        # their signal in the central pixel (peak/aperture ≈ 1.0); real
        # PSFs spread it (peak/aperture ≈ 0.2-0.4 for FWHM≈3 and
        # aperture_radius=2*FWHM).
        if self.fs_max_concentration is not None and n_cand > 0:
            denom = np.where(aperture_excess > 0, aperture_excess, 1.0)
            concentration = peak_excess / denom
            keep = (aperture_excess > 0) & (concentration <= self.fs_max_concentration)
            coo = coo[keep]
            val = val[keep]
            aperture_excess = aperture_excess[keep]
            peak_excess = peak_excess[keep]

        # Pixel-count cull — count mask1 pixels inside each candidate's
        # aperture box. A PSF spans many; an isolated noise spike spans 1.
        if min_pixels_above_threshold is not None and len(coo) > 0:
            min_px = max(1, int(round(int(min_pixels_above_threshold) / b ** 2)))   # piksele zbinowane
            counts = np.zeros(len(coo), dtype=int)
            for i in range(len(coo)):
                yi, xi = int(coo[i, 0]), int(coo[i, 1])
                y0 = max(0, yi - ap_radius); y1 = min(H, yi + ap_radius + 1)
                x0 = max(0, xi - ap_radius); x1 = min(W, xi + ap_radius + 1)
                counts[i] = int(mask1[y0:y1, x0:x1].sum())
            keep = counts >= min_px
            coo = coo[keep]
            val = val[keep]
            aperture_excess = aperture_excess[keep]
            peak_excess = peak_excess[keep]

        self.stats["stars"] = {}
        if len(coo) > 0:
            if self.fs_rank_by == "smoothed":
                sort_keys = data2[coo[:, 0], coo[:, 1]]
            elif self.fs_rank_by == "aperture":
                sort_keys = aperture_excess
            else:
                sort_keys = val
            sorted_i = np.argsort(sort_keys)[::-1]
            self.coo = coo[sorted_i]
            self.adu = val[sorted_i]
            if b > 1:
                self.coo, self.adu = self._refine_on_original(self.coo, exclude)

            self.stats["stars"] = {
                "x": self.coo[:, 1],
                "y": self.coo[:, 0],
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

        self.stars = Table(self.stats["stars"])

    def calc_frame_fwhm(self,threshold=10, fwhm=10, box=10, N_stars=20, clip=4):
        """Full frame analysis: masks (saturation, cosmics, threshold), sky, lines, stars,
        FWHM and shape -> self.stats["frame"] (fwhm, ellipticity, used_stars, ...).

        Runs mk_saturation_mask, mk_stats, sky_map, mk_cosmics_mask, mk_threshold_mask,
        find_lines, classify_lines, find_stars, star_info and calc_star_stats with default
        parameters except the ones below. used_stars tells how many stars the statistics
        come from (0 -> NaN).

        Binning is set in the constructor. Without it (default) large frames, e.g.
        4096 x 4096, may run out of memory in find_lines() - use FFS(image, autobin=True).

        Parameters
        ----------
        threshold : float
            Star detection threshold in units of background sigma (find_stars).
        fwhm : float
            Adopted PSF FWHM in pixels (find_stars).
        box : int
            Half-size of the measurement box around each star, in pixels (star_info).
        N_stars : int | None
            Number of the brightest good stars to measure; None = all (star_info).
        clip : float
            Sigma clipping of outlying stars in the frame statistics; 0 = none (calc_star_stats).
        """
        self.mk_saturation_mask()
        self.mk_stats()
        self.sky_map()
        self.mk_cosmics_mask()
        self.mk_threshold_mask()
        self.find_lines()
        self.classify_lines()
        self.find_stars(threshold=threshold, fwhm=fwhm)
        self.star_info(box=box,N_stars=N_stars)
        self.calc_star_stats(clip=clip)

    def _clip_mask(self,x,clip=4):
        med = np.nanmedian(x)
        sig = mad_std(x, ignore_nan=True)

        if not np.isfinite(sig) or sig == 0:
            return np.ones_like(x, dtype=bool)

        return np.abs(x - med) < clip * sig

    @staticmethod
    def _axial_spread(theta):
        # Circular spread of orientations (period pi): theta and theta +/- pi
        # are the same axis, so a plain std would explode on sign flips.
        # Equals the ordinary std (radians) for small spreads.
        theta = np.asarray(theta, dtype=float)
        theta = theta[np.isfinite(theta)]
        if theta.size == 0:
            return np.nan
        R = np.abs(np.mean(np.exp(2j * theta)))
        return 0.5 * np.sqrt(-2.0 * np.log(max(R, 1e-12)))

    def calc_star_stats(self,clip=4):
        """Frame PSF statistics from the stars measured by star_info() -> self.stats["frame"]
        (fwhm, fwhm_x, fwhm_y, ellipticity, theta_spread, cpe, shape, ci, used_stars).

        Medians over stars; stars whose log(FWHM), log(CPE) or ellipticity deviate by more
        than clip * mad_std are rejected first (if fewer than 3 stars remain, all valid stars are used).

        Parameters
        ----------
        clip : float
            Sigma clipping level in units of mad_std; 0 = no clipping.
        """

        mk = np.isfinite(self.fwhm)

        if clip:
            fwhm = np.where(self.fwhm > 0, np.log(self.fwhm), np.nan)
            cpe = np.where(self.cpe > 0, np.log(self.cpe), np.nan)
            ell = self.ellipticity

            valid = np.isfinite(fwhm) & np.isfinite(cpe) & np.isfinite(ell)

            mk = (valid & self._clip_mask(fwhm, clip) & self._clip_mask(cpe, clip) & self._clip_mask(ell, clip))

            if np.sum(mk) < 3:
                mk = valid

        self.frame_fwhm = np.nanmedian(self.fwhm[mk])
        self.frame_fwhm_x = np.nanmedian(self.fwhm_x[mk])
        self.frame_fwhm_y = np.nanmedian(self.fwhm_y[mk])
        self.frame_ellipticity = np.nanmedian(self.ellipticity[mk])
        self.frame_theta_spread = FFS._axial_spread(self.theta[mk])
        self.frame_cpe = np.nanmedian(self.cpe[mk])
        self.frame_shape = np.nanmedian(self.shape[mk])
        self.frame_ci = np.nanmedian(self.ci[mk])
        self.frame_used_stars = len(self.fwhm[mk])


        self.stats["frame"].update({
            "fwhm": self.frame_fwhm,
            "fwhm_x": self.frame_fwhm_x,
            "fwhm_y": self.frame_fwhm_y,
            "ellipticity": self.frame_ellipticity,
            "theta_spread": self.frame_theta_spread,
            "cpe": self.frame_cpe,
            "shape": self.frame_shape,
            "ci": self.frame_ci,
            "used_stars": self.frame_used_stars
        })


        self.stats_description["frame"].update({
            "fwhm": "Median of the full width at half maximum (FWHM) of the brightest sources in the frame",
            "fwhm_x": "Median FWHM measured along the X axis",
            "fwhm_y": "Median FWHM measured along the Y axis",
            "ellipticity": "Median source ellipticity (1 − b/a), describing PSF elongation",
            "theta_spread": "Axial (period pi) circular spread of source position angles, in radians. Low spread with high ellipticity indicates a common elongation direction (e.g. tracking or wind)",
            "cpe": "Median central pixel excess (CPE) of detected sources; higher values indicate sharper, more centrally concentrated profiles",
            "shape": "Median PSF shape parameter defined as CPE × FWHM^2 ",
            "ci": "Median concentration index (CI) of detected stars in the frame; ",
            "used_stars": "How many star were used to evaluate thos statistics"
        })

    def mk_star_maps(self, scale=0.2, grid=10, min_stars=5, nodes=65):
        """Star maps: PSF from the stars measured by star_info() -> maps["psf_fwhm"],
        maps["psf_ellipticity"], maps["psf_theta"], maps["psf_weight"]; star density from all
        detections in self.coo -> maps["star_density"]; and self.psf_grid (table for a whisker plot).

        No model: at every point a Gaussian-weighted mean of the stars, weight exp(-d^2 / 2 s^2),
        s = scale * longer frame side. theta has period pi, so the components e1 = e cos 2theta,
        e2 = e sin 2theta are averaged and ellipticity = hypot(e1, e2), theta = atan2(e2, e1) / 2
        (radians, from the x axis). Random orientations average out, so the ellipticity map shows
        the coherent part of the elongation. psf_weight is the sum of the weights (about the number
        of stars within ~s); where it is small the map is extrapolated from far stars.
        star_density: Gaussian kernel density of all detections (also saturated and not measured),
        same s, in stars per 100 x 100 px; near the frame edges divided by the part of the kernel
        inside the frame (otherwise the density would drop at the edges). Depends on the detection
        threshold of find_stars().
        The means are computed on nodes x nodes points and linearly interpolated to the image size.

        Parameters
        ----------
        scale : float
            Kernel sigma as a fraction of the longer frame side.
        grid : int
            The whisker table self.psf_grid has grid x grid points (centres of grid x grid cells).
        min_stars : int
            Minimal number of measured stars; fewer -> ValueError.
        nodes : int
            Number of nodes along each axis of the coarse grid.
        """
        ok = np.isfinite(self.fwhm) & np.isfinite(self.ellipticity) & np.isfinite(self.theta)
        if ok.sum() < min_stars:
            raise ValueError(f"mk_star_maps: {ok.sum()} measured stars, at least {min_stars} needed (run star_info)")
        py, px = self.coo[ok, 0].astype(float), self.coo[ok, 1].astype(float)
        e, t = self.ellipticity[ok], self.theta[ok]
        values = np.column_stack((self.fwhm[ok], e * np.cos(2 * t), e * np.sin(2 * t)))

        ny, nx = self.image.shape
        s = scale * max(nx, ny)

        # wezly rzadkiej siatki (z krawedziami kadru) -> pelny rozmiar interpolacja liniowa;
        # theta liczona dopiero z e1, e2 w pelnym rozmiarze (kat ma okres pi)
        gx, gy = np.meshgrid(np.linspace(0, nx - 1, nodes), np.linspace(0, ny - 1, nodes))
        means, weight = FFS._kernel_smooth(px, py, values, gx, gy, s)
        density = self._star_density(gx, gy, s)
        yy, xx = np.mgrid[0:ny, 0:nx]
        coords = np.array([yy * (nodes - 1) / max(ny - 1, 1), xx * (nodes - 1) / max(nx - 1, 1)])
        fwhm, e1, e2, w, dens = (map_coordinates(a, coords, order=1)
                                 for a in (means[..., 0], means[..., 1], means[..., 2], weight, density))

        self.maps["psf_fwhm"] = fwhm
        self.maps["psf_ellipticity"] = np.hypot(e1, e2)
        self.maps["psf_theta"] = 0.5 * np.arctan2(e2, e1)
        self.maps["psf_weight"] = w
        self.maps["star_density"] = dens
        kernel = f"Gaussian-weighted mean of {ok.sum()} stars, sigma {s:.0f} px (mk_star_maps)"
        self.maps_description["psf_fwhm"] = f"PSF FWHM (px): {kernel}"
        self.maps_description["psf_ellipticity"] = f"PSF ellipticity hypot(e1, e2), coherent part: {kernel}"
        self.maps_description["psf_theta"] = f"PSF position angle atan2(e2, e1) / 2 (rad, from x axis): {kernel}"
        self.maps_description["psf_weight"] = f"Sum of kernel weights (about the number of stars within sigma): {kernel}"
        self.maps_description["star_density"] = (
            f"Star density (detections per 100 x 100 px): Gaussian kernel density of {len(self.coo)} detections, "
            f"sigma {s:.0f} px, edge corrected (mk_star_maps)"
        )

        cx = (np.arange(grid) + 0.5) * nx / grid
        cy = (np.arange(grid) + 0.5) * ny / grid
        gx, gy = np.meshgrid(cx, cy)
        means, weight = FFS._kernel_smooth(px, py, values, gx, gy, s)
        self.psf_grid = Table({
            "x": gx.ravel(), "y": gy.ravel(),
            "fwhm": means[..., 0].ravel(),
            "ellipticity": np.hypot(means[..., 1], means[..., 2]).ravel(),
            "theta": 0.5 * np.arctan2(means[..., 2], means[..., 1]).ravel(),
            "weight": weight.ravel(),
            "density": self._star_density(gx, gy, s).ravel(),
        })

    def _star_density(self, gx, gy, s):
        """Gaussian kernel density of all detections (self.coo) at points gx, gy, stars per 100 x 100 px,
        divided by the part of the kernel inside the frame."""
        ny, nx = self.image.shape
        coo = np.asarray(self.coo, dtype=float).reshape(-1, 2)
        if len(coo) == 0:
            return np.zeros(gx.shape)
        ones = np.ones((len(coo), 1))
        _, W = FFS._kernel_smooth(coo[:, 1], coo[:, 0], ones, gx, gy, s)
        inside_x = 0.5 * (erf(gx / (np.sqrt(2) * s)) + erf((nx - 1 - gx) / (np.sqrt(2) * s)))
        inside_y = 0.5 * (erf(gy / (np.sqrt(2) * s)) + erf((ny - 1 - gy) / (np.sqrt(2) * s)))
        return W / (2 * np.pi * s ** 2) / (inside_x * inside_y) * 1e4

    def calc_star_field_stats(self, min_stars=5):
        """Global PSF field statistics from the stars measured by star_info() -> stats["frame"].

        Ellipticity components of each star: e1 = e cos 2theta, e2 = e sin 2theta. Relative to
        the frame centre (phi - position angle of the star seen from the centre):
        - psf_e1, psf_e2, psf_e_mean, psf_theta_mean: mean e1, e2 and the common elongation
          (all stars elongated in one direction: tracking, wind),
        - psf_e_t = mean of -(e1 cos 2phi + e2 sin 2phi): > 0 tangential pattern (field rotation,
          astigmatism), < 0 radial pattern,
        - psf_e_x = mean of (e1 sin 2phi - e2 cos 2phi): pattern rotated by 45 deg (swirl).
        FWHM fitted with a + b x + c y + d r^2 (x, y, r from the frame centre):
        - psf_fwhm_tilt = hypot(b, c) * longer frame side: FWHM change across the frame
          (px per frame, e.g. detector tilt), psf_fwhm_tilt_angle = atan2(c, b): direction of
          growing FWHM (rad, from the x axis),
        - psf_fwhm_curv = d * (half diagonal)^2: FWHM change from the centre to a corner
          (px, e.g. field curvature).
        The optical axis is assumed at the frame centre.

        Parameters
        ----------
        min_stars : int
            Minimal number of measured stars; fewer -> ValueError.
        """
        ok = np.isfinite(self.fwhm) & np.isfinite(self.ellipticity) & np.isfinite(self.theta)
        if ok.sum() < min_stars:
            raise ValueError(f"calc_star_field_stats: {ok.sum()} measured stars, at least {min_stars} needed (run star_info)")
        ny, nx = self.image.shape
        x = self.coo[ok, 1] - (nx - 1) / 2
        y = self.coo[ok, 0] - (ny - 1) / 2
        e, t, fwhm = self.ellipticity[ok], self.theta[ok], self.fwhm[ok]

        e1, e2 = e * np.cos(2 * t), e * np.sin(2 * t)
        phi = np.arctan2(y, x)
        c2, s2 = np.cos(2 * phi), np.sin(2 * phi)
        m1, m2 = np.mean(e1), np.mean(e2)

        r2 = x ** 2 + y ** 2
        A = np.column_stack((np.ones_like(x), x, y, r2))
        (_, b, c, d), *_ = np.linalg.lstsq(A, fwhm, rcond=None)

        psf = {
            "psf_e1": m1,
            "psf_e2": m2,
            "psf_e_mean": np.hypot(m1, m2),
            "psf_theta_mean": 0.5 * np.arctan2(m2, m1),
            "psf_e_t": np.mean(-(e1 * c2 + e2 * s2)),
            "psf_e_x": np.mean(e1 * s2 - e2 * c2),
            "psf_fwhm_tilt": np.hypot(b, c) * max(nx, ny),
            "psf_fwhm_tilt_angle": np.arctan2(c, b),
            "psf_fwhm_curv": d * ((nx - 1) ** 2 + (ny - 1) ** 2) / 4,
        }
        self.stats["frame"].update(psf)
        self.stats_description["frame"].update({
            "psf_e1": "Mean e1 = e cos 2theta of the measured stars",
            "psf_e2": "Mean e2 = e sin 2theta of the measured stars",
            "psf_e_mean": "Common elongation hypot(mean e1, mean e2) (tracking, wind)",
            "psf_theta_mean": "Direction of the common elongation (rad, from the x axis)",
            "psf_e_t": "Mean tangential ellipticity around the frame centre (> 0 tangential, < 0 radial)",
            "psf_e_x": "Mean cross ellipticity around the frame centre (45 deg swirl)",
            "psf_fwhm_tilt": "FWHM change across the frame from a linear + r^2 fit (px per frame)",
            "psf_fwhm_tilt_angle": "Direction of growing FWHM (rad, from the x axis)",
            "psf_fwhm_curv": "FWHM change from the frame centre to a corner, r^2 term (px)",
        })

    @staticmethod
    def _kernel_smooth(px, py, values, gx, gy, s):
        """Gaussian-weighted means of values (n_stars x n_columns) at points gx, gy (any shape)
        -> (means of shape gx.shape + (n_columns,), sum of weights of shape gx.shape)."""
        d2 = (gx.reshape(-1, 1) - px) ** 2 + (gy.reshape(-1, 1) - py) ** 2
        w = np.exp(-0.5 * d2 / s ** 2)
        W = w.sum(axis=1)
        means = (w @ values) / W[:, None]
        return means.reshape(gx.shape + (values.shape[1],)), W.reshape(gx.shape)

    def star_info(self, box=10, N_stars=None, exclude=EXCLUDE_MEASURE):
        """Measure the stars found by find_stars() on the original image (no binning): a box of
        +-box px around each position, FWHM, ellipticity, theta, CPE, CI -> self.stars / stats["stars"].

        Measures the first (brightest) N_stars stars, skipping saturated ones and those with
        bad_mask() pixels in the box. All detections stay in self.coo / self.stars; the stars
        that were not measured have NaN measurements. Warns when the median FWHM > 0.8 * box (the box cuts the profile).

        Parameters
        ----------
        box : int
            Half-size of the measurement box, in pixels. Should be >= 1.5 * FWHM.
        N_stars : int | None
            Number of stars to measure; None = all detected stars.
        exclude : sequence of str
            Names of the masks whose pixels are skipped (bad_mask); missing masks are ignored.
            Without "stars" by default, so that the measured stars are not masked.
        """

        n = len(self.coo)

        self.box_mag = np.full(n, np.nan)
        self.bkg = np.full(n, np.nan)
        self.ellipticity = np.full(n, np.nan)
        self.theta = np.full(n, np.nan)
        self.fwhm = np.full(n, np.nan)
        self.fwhm_x = np.full(n, np.nan)
        self.fwhm_y = np.full(n, np.nan)
        self.cpe = np.full(n, np.nan)
        self.shape = np.full(n, np.nan)
        self.ci = np.full(n, np.nan)

        if N_stars is None:
            N_stars = n

        # find_stars bez detekcji zostawia listy; mierzymy nawet jedna gwiazde (used_stars mowi ile)
        self.coo = np.asarray(self.coo, dtype=int).reshape(-1, 2)
        self.adu = np.asarray(self.adu, dtype=float)

        # wszystkie detekcje zostaja w coo (np. do maski gwiazd); niezmierzone maja NaN
        bad = self.bad_mask(exclude)
        n_measured = 0

        for i, (y, x) in enumerate(self.coo):

            if n_measured >= N_stars:
                break

            if self.adu[i] >= self.saturation:
                continue

            y0 = max(0, y - box)
            y1 = min(self.image.shape[0], y + box)
            x0 = max(0, x - box)
            x1 = min(self.image.shape[1], x + box)

            cut = self.image[y0:y1, x0:x1]

            if cut.size == 0:
                continue

            if bad[y0:y1, x0:x1].any():  # zle piksele w wycinku psuja pomiar (fotometria: photutils)
                continue

            top = cut[:1, :]
            bottom = cut[-1:, :]
            left = cut[:, :1]
            right = cut[:, -1:]

            bkg_pixels = np.concatenate([top.ravel(),bottom.ravel(),left.ravel(),right.ravel()])

            bkg = np.median(bkg_pixels)
            cut = cut - bkg

            flux = np.sum(cut)
            if flux <= 0:
                continue

            box_mag = -2.5 * np.log10(flux) + 25

            # Adaptive moments resist noise and neighbours in the box; plain
            # pca() is the fallback when they do not converge (e.g. donuts).
            _, e, t = FFS.adaptive_moments(cut)
            if not np.isfinite(e):
                _, e, t = FFS.pca(cut)

            fx, fy = FFS.fwhm(cut)
            if np.isnan(fx) or np.isnan(fy):
                fx = fy = star_fwhm = np.nan
            else:
                star_fwhm = (fx + fy) / 2

            cpe = FFS.cpe(cut)
            shape = cpe * star_fwhm**2

            # Radii scaled to this star's own measured FWHM; fall back to
            # fixed radii when the FWHM measurement failed (NaN / non-positive).
            if np.isfinite(star_fwhm) and star_fwhm > 0:
                r1 = 1.0 * star_fwhm
                r2 = 2.0 * star_fwhm
            else:
                r1 = 2.0
                r2 = 4.0

            ci = FFS.concentration_index(cut,r1,r2)

            if not (np.isfinite(star_fwhm) and np.isfinite(e) and np.isfinite(ci)):
                continue

            n_measured += 1
            (self.box_mag[i], self.bkg[i], self.ellipticity[i], self.theta[i], self.fwhm[i],
             self.fwhm_x[i], self.fwhm_y[i], self.cpe[i], self.shape[i], self.ci[i]) = (
                box_mag, bkg, e, t, star_fwhm, fx, fy, cpe, shape, ci)

        self.stats["stars"]["box_mag"] = self.box_mag
        self.stats["stars"]["bkg"] = self.bkg
        self.stats["stars"]["fwhm"] = self.fwhm
        self.stats["stars"]["fwhm_x"] = self.fwhm_x
        self.stats["stars"]["fwhm_y"] = self.fwhm_y
        self.stats["stars"]["ellipticity"] = self.ellipticity
        self.stats["stars"]["theta"] = self.theta
        self.stats["stars"]["cpe"] = self.cpe
        self.stats["stars"]["shape"] = self.shape
        self.stats["stars"]["ci"] = self.ci
        self.stats["stars"]["x"] = self.coo[:, 1]
        self.stats["stars"]["y"] = self.coo[:, 0]
        self.stats["stars"]["max_adu"] = self.adu


        self.stars = Table(self.stats["stars"])

        # wycinek ma +-box px; przy FWHM porownywalnym z box obcina skrzydla i zaniza FWHM
        # (synth.: -4% przy FWHM = box, < 1% przy FWHM = 0.8 box)
        if len(self.fwhm) and np.nanmedian(self.fwhm) > 0.8 * box:
            warnings.warn(
                f"FFS.star_info: median FWHM {np.nanmedian(self.fwhm):.1f} px is comparable to box={box}; "
                f"FWHM may be underestimated, use a larger box (e.g. box >= {int(np.ceil(1.5 * np.nanmedian(self.fwhm)))})"
            )

        self.stats_description["stars"].update({

            "box_mag": (
                "magnitude of star computed as the sum in the box - background "
            ),

            "bkg": (
                "Local sky background value estimated from the border pixels "
                "of the analysis box (in ADU); calculated as the median of "
                "all edge pixels to minimize stellar flux contamination"
            ),

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
                "Ellipticity of each detected star, "
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

            "shape": (
                "PSF shape parameter defined as CPE × FWHM^2; "
                "for a Gaussian PSF this is approximately constant, "
                "deviations indicate non-Gaussian profiles (e.g. broader wings or sharper cores)"
            ),

            "ci": (
                "Concentration Index (CI), defined as the ratio of flux within inner radius r1 "
                "to flux within outer radius r2 (CI = F(r<r1)/F(r<r2)); "
                "higher values indicate more centrally concentrated profiles"
            ),

        })

    def sky_map(self,n_segments=10, exclude=EXCLUDE_ALL):
        """Sky background map -> maps["sky"], original size.

        Medians of n_segments x n_segments segments (on the binned image, without bad_mask())
        linearly interpolated between segment centres. The band between the outermost centres
        and the frame edges is linearly extrapolated: points on the frame edges get values
        extrapolated from the two nearest segments of their row / column (corners from the edge
        points), then all points are interpolated together. No model (polynomial) in the map.
        Calls sky_gradient(), so also fills its statistics.

        Parameters
        ----------
        n_segments : int
            Number of segments along each axis.
        exclude : sequence of str
            Names of the masks whose pixels are skipped (bad_mask); missing masks are ignored.
        """
        self.sky_gradient(n_segments=n_segments, exclude=exclude)

        ny, nx = self.image.shape

        # siatka srodkow segmentow (NaN - segment pominiety) obramowana krawedziami kadru:
        # najpierw kolumny x = 0 i nx-1 (ekstrapolacja wzdluz wierszy), potem wiersze y = 0
        # i ny-1 (wzdluz kolumn, razem z nowymi kolumnami, co daje narozniki)
        sx = np.asarray(self.stats["sky"]["surface_x"])
        sy = np.asarray(self.stats["sky"]["surface_y"])
        xs = np.unique(sx)
        ys = np.unique(sy)
        X = np.concatenate(([0], xs, [nx - 1]))
        Y = np.concatenate(([0], ys, [ny - 1]))
        grid = np.full((len(Y), len(X)), np.nan)
        grid[np.searchsorted(ys, sy) + 1, np.searchsorted(xs, sx) + 1] = self.stats["sky"]["surface_val"]
        for i in range(1, len(Y) - 1):
            grid[i, 0] = FFS._extrapolate(X[1:-1], grid[i, 1:-1], X[0])
            grid[i, -1] = FFS._extrapolate(X[1:-1], grid[i, 1:-1], X[-1])
        for j in range(len(X)):
            grid[0, j] = FFS._extrapolate(Y[1:-1], grid[1:-1, j], Y[0])
            grid[-1, j] = FFS._extrapolate(Y[1:-1], grid[1:-1, j], Y[-1])

        gx, gy = np.meshgrid(X, Y)
        ok = np.isfinite(grid)
        points = np.column_stack((gx[ok], gy[ok]))
        values = grid[ok]

        yy, xx = np.mgrid[0:ny, 0:nx]

        sky_map = LinearNDInterpolator(points, values)(xx, yy)

        # poza otoczka wypukla punktow (tylko gdy brakuje skrajnych segmentow) LinearNDInterpolator daje NaN
        outside = np.isnan(sky_map)
        if outside.any():
            sky_map[outside] = NearestNDInterpolator(points, values)(xx[outside], yy[outside])

        self.maps["sky"] = sky_map
        self.maps_description["sky"] = (
            "Sky background map: linear interpolation (Delaunay triangulation) of segment medians, "
            "linearly extrapolated to the frame edges"
        )

    def sky_gradient(self,n_segments=10, min_good_fraction=0.5, exclude=EXCLUDE_ALL):
        """Sky gradient: medians of n_segments x n_segments segments fitted with a 2nd order
        polynomial surface -> stats["frame"] (bkg_max_amplitude, bkg_frame_gradient,
        sky_surface_coeff), stats["sky"] and self.sky (table of segment centres and fitted values).
        bkg_max_amplitude: max - min of the model over the frame (ADU); bkg_frame_gradient: max
        slope of the model over the frame times the longer frame side (ADU per frame).

        Parameters
        ----------
        n_segments : int
            Number of segments along each axis.
        min_good_fraction : float
            Segments with a smaller fraction of good pixels (outside bad_mask()) are skipped.
        exclude : sequence of str
            Names of the masks whose pixels are skipped (bad_mask); missing masks are ignored.
        """
        # mediany segmentow liczone na obrazie zbinowanym; srodki segmentow (i wspolczynniki
        # wielomianu) przeliczone na wspolrzedne oryginalu
        image = self.binned_image
        segments = FFS.make_segments(image, n_segments=n_segments)
        bad_segments = FFS.make_segments(self.to_binned_mask(self.bad_mask(exclude)), n_segments=n_segments)

        x = []
        y = []
        back = []
        back_mad = []
        for s, b in zip(segments, bad_segments):
            good = s["subframe"][~b["subframe"]]
            if good.size < min_good_fraction * s["subframe"].size:
                continue                   # za malo dobrych pikseli (bad_mask) w segmencie
            x_tmp = (s["x"][0] + s["x"][1]) // 2
            y_tmp = (s["y"][0] + s["y"][1]) // 2
            x.append(x_tmp)
            y.append(y_tmp)
            back.append(np.median(good))
            back_mad.append(mad_std(good))

        x = np.rint(self.to_original_coord(x)).astype(int)   # calkowite, jak przed binowaniem
        y = np.rint(self.to_original_coord(y)).astype(int)
        back = np.array(back)

        A = np.vstack([np.ones_like(x), x, y, x ** 2, x * y, y ** 2]).T
        coeff, *_ = np.linalg.lstsq(A, back, rcond=None)

        bk = [FFS.polysurf(xi, yi, coeff) for xi, yi in zip(x, y)]

        # parametry z modelu na calym kadrze (oryginal): amplituda = max - min na siatce
        # 101 x 101 (lapie tez ekstremum wewnatrz, np. srodek winietowania); gradient
        # wielomianu 2. stopnia jest liniowy, wiec max |grad| jest w ktoryms narozniku
        ny, nx = self.image.shape
        gx, gy = np.meshgrid(np.linspace(0, nx - 1, 101), np.linspace(0, ny - 1, 101))
        surf = FFS.polysurf(gx, gy, coeff)
        max_amplitude = surf.max() - surf.min()

        a0, a1, a2, a3, a4, a5 = coeff
        corners = np.array([[0, 0], [nx - 1, 0], [0, ny - 1], [nx - 1, ny - 1]], dtype=float)
        dx = a1 + 2 * a3 * corners[:, 0] + a4 * corners[:, 1]
        dy = a2 + a4 * corners[:, 0] + 2 * a5 * corners[:, 1]
        frame_gradient = np.hypot(dx, dy).max() * max(nx, ny)   # ADU na kadr

        self.max_amplitude = max_amplitude
        self.frame_gradient = frame_gradient
        self.sky_surface_coeff = coeff
        self.sky_surface_bkg = bk
        self.sky_surface_x = x
        self.sky_surface_y = y

        self.stats["frame"].update({
            "bkg_max_amplitude": self.max_amplitude,
            "bkg_frame_gradient": self.frame_gradient,
            "sky_surface_coeff": self.sky_surface_coeff,
        })

        self.stats["sky"] = {
            "surface_bkg": self.sky_surface_bkg,
            "surface_x": self.sky_surface_x,
            "surface_y": self.sky_surface_y,
            "surface_val": back,
        }

        self.sky = Table([self.sky_surface_x,self.sky_surface_y,self.sky_surface_bkg], names=["sky_surface_x","sky_surface_y","sky_surface_bkg"])


        self.stats_description["frame"].update({

            "bkg_max_amplitude": (
                "Max - min of the sky model (2nd order polynomial) over the frame (ADU)"
            ),

            "bkg_frame_gradient": (
                "Max slope |grad| of the sky model over the frame times the longer frame side (ADU per frame)"
            ),

            "sky_surface_coeff": (
                "Coefficients of the fitted sky model describing "
                "the sky background variation across the image"
            ),
        })

        self.stats_description["sky"] = {
            "surface_bkg": (
                "Estimated sky background surface evaluated over the image grid, "
            ),
            "surface_x": (
                "X-coordinate grid used for evaluating the sky background surface "
                "(pixel indices, 0-based)"
            ),
            "surface_y": (
                "Y-coordinate grid used for evaluating the sky background surface "
                "(pixel indices, 0-based)"
            ),
            "surface_val": (
                "Median of each segment at (surface_x, surface_y), in ADU"
            ),
        }



    def find_lines(self, min_length=30, fwhm=4.0, max_gap=5, min_fill=0.5, half_width=3.0,
                   max_candidates=100, steps=None, max_shift=1.0, exclude=EXCLUDE_ALL,
                   ring_radius=8, ring_omax=0.4):
        """Detect trails (lines) with the Hough transform of binned_masks["threshold"]
        -> self.lines (table: val, rho, theta, x0, y0, x1, y1), stats["lines"], masks["lines"].

        Only line-like pixels vote (ring test, see hough_transform); segments are checked on the
        full threshold mask. Candidates are the accumulator cells with most votes. A candidate is
        a line when the mask pixels along it form a continuous segment: gaps <= max_gap, length >= 10 * fwhm,
        fill >= min_fill. Pixels of an accepted line (within half_width) and their votes are
        removed before the next candidate. Interrupted collinear pieces of one trail become one line.
        Lines follow x*cos(theta) + y*sin(theta) = rho (x = column, y = row); (x0, y0)-(x1, y1)
        are the segment ends. Detection runs on the binned image; results are in original pixels.
        Requires mk_threshold_mask() (computed with default nsigma if missing).

        Parameters
        ----------
        min_length : float
            Minimal number of votes (mask pixels on the line) for a candidate, in original pixels.
        fwhm : float
            Trail width scale in original pixels; the minimal segment length is 10 * fwhm.
        max_gap : float
            Longest allowed gap in a segment, in original pixels.
        min_fill : float
            Minimal fraction of the segment covered by mask pixels (0-1).
        half_width : float
            Half-width of the band removed around an accepted line (binned pixels) and of the
            line in masks["lines"] (original pixels).
        max_candidates : int
            Maximal number of accumulator maxima checked.
        steps : int | None
            Number of theta angles in the Hough transform; None chooses it from max_shift.
        max_shift : float
            Maximal shift of the ends of the longest line caused by the theta step, in binned
            pixels (used when steps is None).
        exclude : sequence of str
            Names of the masks whose pixels are unknown (bad_mask), not gaps, when segments
            are checked; missing masks are ignored.
        ring_radius : float | None
            Radius of the ring test in binned pixels; None = all threshold pixels vote (hough_transform).
        ring_omax : float
            Maximal fraction of the ring in the mask (hough_transform).
        """
        self.hough_transform(steps=steps, max_shift=max_shift, ring_radius=ring_radius, ring_omax=ring_omax)

        # dlugosci (min_length, fwhm, max_gap) przeliczone na piksele obrazu zbinowanego;
        # half_width (pas usuwania przyjetej linii, tolerancja jak max_shift) zostaje w pikselach
        # zbinowanych - po podzieleniu przez bin konce dlugich sladow wychodzily z pasa
        b = self.bin
        min_votes = min_length / b
        min_segment = 10 * fwhm / b
        max_gap_b = max(1, int(round(max_gap / b)))
        half_width_b = half_width

        rh, theta = self.rh, self.hough_theta
        rh0 = -int(rh[0])                     # koszyk k <-> rho = k - rh0
        cos_t, sin_t = np.cos(theta), np.sin(theta)
        # glosujace piksele (do odejmowania glosow) i wszystkie piksele maski (do usuwania z maski roboczej)
        vy, vx = (a.astype(float) for a in np.nonzero(self.binned_masks["line_votes"]))
        ys, xs = (a.astype(float) for a in np.nonzero(self.binned_masks["threshold"]))

        # znajdujemy maksimum, jak jest powyzej th, a nastepnie usuwamy z glosowania piksele ktore sa w odleglosci
        # half_width od lini

        acc = self.accumulator.copy()
        work = self.binned_masks["threshold"].copy()   # maska, z ktorej usuwamy piksele przyjetych linii
        # piksele "nieznane": zle (bad_mask) i usuniete razem z przyjetymi liniami
        removed = self.to_binned_mask(self.bad_mask(exclude))
        lines_rho, lines_theta, lines_val = [], [], []
        lines_x0, lines_y0, lines_x1, lines_y1 = [], [], [], []

        for _ in range(max_candidates):
            ri_max, ti_max = np.unravel_index(np.argmax(acc), acc.shape)
            val = acc[ri_max, ti_max]
            if val < min_votes:
                break
            rho, t = rh[ri_max], theta[ti_max]

            # czy wzdluz linii piksele maski tworza ciagly odcinek?
            pos, occupied = FFS._line_profile(work, rho, t)
            _, unknown = FFS._line_profile(removed, rho, t)   # skrzyzowania z juz znalezionymi liniami
            segments = FFS._line_segments(occupied, max_gap_b, min_segment, min_fill, unknown)
            if not segments:
                # odrzucony: wylaczamy te komorke i jej sasiadow, piksele zostaja
                acc[max(0, ri_max - 2):ri_max + 3, max(0, ti_max - 1):ti_max + 2] = -1
                continue

            # przyjety: zapisujemy linie od poczatku pierwszego do konca ostatniego odcinka
            # (wspolliniowe fragmenty przerwanego sladu to jeden slad)
            start, stop = segments[0][0], segments[-1][1]
            c, s = np.cos(t), np.sin(t)
            p0, p1 = pos[start], pos[stop - 1]
            lines_rho.append(rho)
            lines_theta.append(t)
            lines_val.append(val)
            lines_x0.append(rho * c - p0 * s)
            lines_y0.append(rho * s + p0 * c)
            lines_x1.append(rho * c - p1 * s)
            lines_y1.append(rho * s + p1 * c)

            # usuwamy glosy pikseli tej linii i same piksele (tez z maski roboczej)
            on_votes = np.abs(vx * c + vy * s - rho) <= half_width_b
            FFS._hough_votes(acc, vx[on_votes], vy[on_votes], cos_t, sin_t, rh0, weight=-1)
            vx, vy = vx[~on_votes], vy[~on_votes]
            on_line = np.abs(xs * c + ys * s - rho) <= half_width_b
            work[ys[on_line].astype(int), xs[on_line].astype(int)] = False
            removed[ys[on_line].astype(int), xs[on_line].astype(int)] = True
            xs, ys = xs[~on_line], ys[~on_line]

        # z ukladu zbinowanego do oryginalu: x_o = b*x + h, y_o = b*y + h (h = (b-1)/2), wiec
        # x_o*cos + y_o*sin = b*rho + h*(cos + sin); theta sie nie zmienia
        lines_theta = np.array(lines_theta)
        h = (b - 1) / 2.0
        lines_rho = b * np.array(lines_rho) + h * (np.cos(lines_theta) + np.sin(lines_theta))
        lines_x0, lines_y0 = self.to_original_coord(lines_x0), self.to_original_coord(lines_y0)
        lines_x1, lines_y1 = self.to_original_coord(lines_x1), self.to_original_coord(lines_y1)

        self.lines_rho = np.array(lines_rho)
        self.lines_theta = np.array(lines_theta)
        self.lines_val = np.array(lines_val)
        self.lines_x0 = np.array(lines_x0)
        self.lines_y0 = np.array(lines_y0)
        self.lines_x1 = np.array(lines_x1)
        self.lines_y1 = np.array(lines_y1)

        tmp = {}
        tmp["val"] = self.lines_val
        tmp["rho"] = self.lines_rho
        tmp["theta"] = self.lines_theta
        tmp["x0"] = self.lines_x0
        tmp["y0"] = self.lines_y0
        tmp["x1"] = self.lines_x1
        tmp["y1"] = self.lines_y1
        self.stats["lines"] = tmp

        # tabela linii (do rysowania np. w FitsView) i maska pikseli linii
        self.lines = Table(self.stats["lines"])
        mask = np.zeros(self.image.shape, dtype=bool)
        for row in self.lines:
            mask |= FFS.line_to_mask(mask.shape, row["rho"], row["theta"], half_width,
                              row["x0"], row["y0"], row["x1"], row["y1"])
        self.masks["lines"] = mask
        self.masks_description["lines"] = (
            f"Pixels within {half_width} px of the detected line segments (find_lines)"
        )

        self.stats_description["lines"] = {

            "val": (
                "Detection strength (number of pixels in line), "
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

            "x0": "X of the first end of the detected line segment (pixels, 0-based)",
            "y0": "Y of the first end of the detected line segment (pixels, 0-based)",
            "x1": "X of the second end of the detected line segment (pixels, 0-based)",
            "y1": "Y of the second end of the detected line segment (pixels, 0-based)",
        }

    def classify_lines(self, half_width=3.0, margin=10, min_side=10, min_decay=2.0,
                       exclude=EXCLUDE_LINES):
        """Separate diffraction spikes of bright stars from trails. Requires find_lines().

        Step 1: is there a saturated star (image >= self.saturation) on the segment (with a
        margin at the ends)? If not - a trail.
        Step 2: does the brightness (image - sky) along the line fall with the distance from
        the star? On each side of the star where the segment has >= min_side px, the median
        of the nearer half is divided by the median of the farther half. A spike when the
        ratio is >= min_decay on every side. When no side can be measured, step 1 decides (spike).

        Adds columns spike, star_x, star_y, decay to self.lines (and stats["lines"]), creates
        masks["spikes"] and rebuilds masks["lines"] with trails only.

        Parameters
        ----------
        half_width : float
            Half-width of the lines in masks["spikes"] and masks["lines"], in pixels.
        margin : float
            Margin at the segment ends where a saturated star is still looked for, in pixels.
        min_side : int
            Minimal segment length on one side of the star needed to measure the decay, in pixels.
        min_decay : float
            Minimal near/far brightness ratio for a spike.
        exclude : sequence of str
            Names of the masks whose pixels are skipped (bad_mask); missing masks are ignored.
        """
        level, _, _ = self._background()
        sat = self.image >= self.saturation
        resid = (self.image - level).astype(float)
        resid[sat | self.bad_mask(exclude)] = np.nan

        spike, star_x, star_y, decay = [], [], [], []
        for row in self.lines:
            rho, t = row["rho"], row["theta"]
            c, s = np.cos(t), np.sin(t)
            p0, p1 = sorted([-row["x0"] * s + row["y0"] * c, -row["x1"] * s + row["y1"] * c])

            # etap 1: nasycona gwiazda na odcinku (najdluzszy ciag nasyconych pikseli)
            pos, occupied = FFS._line_profile(sat, rho, t)
            near = occupied & (pos >= p0 - margin) & (pos <= p1 + margin)
            if not near.any():
                spike.append(False)
                star_x.append(np.nan)
                star_y.append(np.nan)
                decay.append(np.nan)
                continue
            labels, _ = label(near)
            idx = np.nonzero(labels == np.argmax(np.bincount(labels)[1:]) + 1)[0]
            pc = pos[idx].mean()
            star_x.append(rho * c - pc * s)
            star_y.append(rho * s + pc * c)

            # etap 2: spadek jasnosci z odlegloscia od gwiazdy, osobno po obu stronach
            vals = FFS._line_values(resid, rho, t, pos)
            ratios = []
            for side in (-1, 1):
                d = side * (pos - pc)
                ok = (d > 0) & (pos >= p0) & (pos <= p1) & np.isfinite(vals)
                if ok.sum() < min_side:
                    continue
                d, v = d[ok], vals[ok]
                half = (d.min() + d.max()) / 2.0
                v_near, v_far = np.median(v[d <= half]), np.median(v[d > half])
                if v_far > 0:
                    ratios.append(v_near / v_far)
                else:                              # dalsza polowa juz na poziomie tla
                    ratios.append(np.inf if v_near > 0 else np.nan)
            spike.append(all(r >= min_decay for r in ratios))   # pusta lista -> etap 1
            decay.append(min(ratios) if ratios else np.nan)

        for name, values in (("spike", spike), ("star_x", star_x), ("star_y", star_y), ("decay", decay)):
            self.stats["lines"][name] = np.array(values)
            self.lines[name] = np.array(values)

        # maski od nowa: promienie do masks["spikes"], w masks["lines"] zostaja tylko slady
        spikes = np.zeros(self.image.shape, dtype=bool)
        lines = np.zeros(self.image.shape, dtype=bool)
        for row in self.lines:
            m = FFS.line_to_mask(spikes.shape, row["rho"], row["theta"], half_width,
                                 row["x0"], row["y0"], row["x1"], row["y1"])
            if row["spike"]:
                spikes |= m
            else:
                lines |= m
        self.masks["spikes"] = spikes
        self.masks_description["spikes"] = (
            f"Pixels within {half_width} px of lines classified as diffraction spikes (classify_lines)"
        )
        self.masks["lines"] = lines
        self.masks_description["lines"] = (
            f"Pixels within {half_width} px of the detected line segments, "
            f"without diffraction spikes (find_lines, classify_lines)"
        )

        self.stats_description["lines"].update({
            "spike": "True if the line is a diffraction spike of a saturated star (classify_lines)",
            "star_x": "X of the saturated star on the line (pixels, 0-based); NaN if none",
            "star_y": "Y of the saturated star on the line (pixels, 0-based); NaN if none",
            "decay": (
                "Smallest (over both sides of the star) ratio of median brightness in the nearer "
                "and the farther half of the segment; NaN if not measured"
            ),
        })

    def hough_transform(self, steps=None, max_shift=1.0, ring_radius=8, ring_omax=0.4):
        """Hough transform of binned_masks["threshold"] (binned image).

        Stores self.accumulator (rho bins x theta angles), self.rh (rho of the bins, step 1 px)
        and self.hough_theta - all in binned image coordinates. The threshold mask is computed
        (mk_threshold_mask() with default nsigma) if missing.

        Only line-like pixels vote (binned_masks["line_votes"], antipodal ring test, _ring_test):
        a mask pixel votes when two opposite points of a circle of radius ring_radius around it
        are in the mask and less than ring_omax of the circle is in the mask. Stars smaller
        than the circle (empty circle) and round blobs of any size (>= half of the circle)
        do not vote; elongated blobs do when their half-length > ring_radius.

        Parameters
        ----------
        steps : int | None
            Number of theta angles; None chooses it so that the longest line (diagonal) does
            not drift at its ends by more than max_shift.
        max_shift : float
            Maximal drift of the line ends, in binned pixels (used when steps is None).
        ring_radius : float | None
            Radius of the ring test in binned pixels (>= 5); None = all threshold pixels vote.
        ring_omax : float
            Maximal fraction of the circle in the mask (lines ~0.1-0.25, round blobs >= 0.5).
        """
        if "threshold" not in self.binned_masks:   # maska liczona sama, jesli jej jeszcze nie ma
            self.mk_threshold_mask()
        if ring_radius is None:
            self.binned_masks["line_votes"] = self.binned_masks["threshold"]
        else:
            self.binned_masks["line_votes"] = FFS._ring_test(self.binned_masks["threshold"], ring_radius, ring_omax)
        ys, xs = np.nonzero(self.binned_masks["line_votes"])  # piksele ktore glosuja
        xs = xs.astype(float)
        ys = ys.astype(float)

        ny, nx = self.binned_image.shape
        rh0 = int((ny ** 2 + nx ** 2) ** 0.5)

        if steps is None:
            # siatka katow gesta na tyle, zeby najdluzsza linia (przekatna) nie odjechala
            # na koncach o wiecej niz max_shift px: przesuniecie ~ dlugosc * krok / 4
            diag = np.hypot(nx, ny)
            steps = int(np.ceil(np.pi * diag / (4 * max_shift)))
            steps += steps % 2      # parzyste, zeby 0 deg bylo w siatce

        theta = np.deg2rad(np.linspace(-90, 90, steps, endpoint=False))   # robimy siatke katow theta
        cos_t = np.cos(theta)
        sin_t = np.sin(theta)
        rh = np.arange(-rh0, rh0 + 1, dtype=float)

        self.accumulator = np.zeros((len(rh), len(theta)))
        FFS._hough_votes(self.accumulator, xs, ys, cos_t, sin_t, rh0)

        self.hough_theta = theta
        self.rh = rh

    @staticmethod
    def _ring_test(mask, radius, omax):
        """Antipodal ring test: mask pixels with two opposite points of a circle of the given
        radius in the mask and less than omax of the circle in the mask (line-like pixels).
        The circle has K >= 2 pi radius / 2 points (every ~2 px, so a trail of width ~2-3 px
        is always hit)."""
        K = int(np.ceil(np.pi * radius))
        K += K % 2
        phi = 2 * np.pi * np.arange(K) / K
        dx = np.rint(radius * np.cos(phi)).astype(int)
        dy = np.rint(radius * np.sin(phi)).astype(int)
        ny, nx = mask.shape
        shifted = []                                   # shifted[k][p] = mask[p + d_k], 0 poza obrazem
        for i, j in zip(dx, dy):
            sh = np.zeros_like(mask)
            sh[max(0, -j):ny - max(0, j), max(0, -i):nx - max(0, i)] = \
                mask[max(0, j):ny - max(0, -j), max(0, i):nx - max(0, -i)]
            shifted.append(sh)
        occupied = np.zeros(mask.shape, dtype=np.int16)
        for sh in shifted:
            occupied += sh
        pair = np.zeros(mask.shape, dtype=bool)
        for k in range(K // 2):                        # para przeciwlegla: k i k + K/2
            pair |= shifted[k] & shifted[k + K // 2]
        return mask & pair & (occupied < omax * K)

    @staticmethod
    def _hough_votes(acc, xs, ys, cos_t, sin_t, rh0, weight=1):
        # dodaje (weight=1) albo odejmuje (weight=-1) glosy pikseli (xs, ys) w akumulatorze:
        # dla kazdego kata theta piksel glosuje na koszyk rho = round(x*cos + y*sin) (+ rh0)
        R, T = acc.shape
        ri = np.round(xs[:, None] * cos_t[None, :] + ys[:, None] * sin_t[None, :] + rh0).astype(int)
        # zliczamy glosy na komorke: komorka (r, t) ma w splaszczonej tablicy numer r * T + t
        flat = ri * T + np.arange(T)
        acc += weight * np.bincount(flat.ravel(), minlength=R * T).reshape(R, T)

    @staticmethod
    def _line_profile(mask, rho, theta):
        """Maska wzdluz linii x*cos(theta) + y*sin(theta) = rho, krokami po 1 px.

        Zwraca (pos, occupied): polozenia wzdluz linii [px] i czy w danym miejscu
        jest piksel maski (ktorykolwiek z 3 pikseli w poprzek linii).
        """
        ny, nx = mask.shape
        c, s = np.cos(theta), np.sin(theta)
        D = int(np.ceil(np.hypot(nx, ny)))
        pos = np.arange(-D, D + 1)
        x = rho * c - pos * s
        y = rho * s + pos * c
        inside = (x > -0.5) & (x < nx - 0.5) & (y > -0.5) & (y < ny - 0.5)
        pos, x, y = pos[inside], x[inside], y[inside]

        occupied = np.zeros(len(pos), dtype=bool)
        for d in (-1, 0, 1):                       # 3 piksele w poprzek linii
            xi = np.rint(x + d * c).astype(int)
            yi = np.rint(y + d * s).astype(int)
            ok = (xi >= 0) & (xi < nx) & (yi >= 0) & (yi < ny)
            occupied[ok] |= mask[yi[ok], xi[ok]]
        return pos, occupied

    @staticmethod
    def _line_values(image, rho, theta, pos):
        """Wartosci obrazu w polozeniach pos wzdluz linii: mediana 3 pikseli w poprzek linii
        (NaN w obrazie pomijane; NaN gdy wszystkie 3 to NaN lub poza kadrem)."""
        ny, nx = image.shape
        c, s = np.cos(theta), np.sin(theta)
        x = rho * c - pos * s
        y = rho * s + pos * c
        vals = np.full((3, len(pos)), np.nan)
        for k, d in enumerate((-1, 0, 1)):         # 3 piksele w poprzek linii
            xi = np.rint(x + d * c).astype(int)
            yi = np.rint(y + d * s).astype(int)
            ok = (xi >= 0) & (xi < nx) & (yi >= 0) & (yi < ny)
            vals[k, ok] = image[yi[ok], xi[ok]]
        with warnings.catch_warnings():
            warnings.simplefilter("ignore", RuntimeWarning)   # All-NaN slice
            return np.nanmedian(vals, axis=0)

    @staticmethod
    def _line_segments(occupied, max_gap=5, min_length=40, min_fill=0.5, unknown=None):
        """Ciagle odcinki wzdluz linii: przerwy <= max_gap, dlugosc >= min_length,
        zajetosc >= min_fill. Zwraca liste (start, stop) - indeksy w occupied.

        unknown: miejsca, z ktorych usunieto piksele juz znalezionych linii
        (skrzyzowania). Nie przerywaja odcinka, ale nie licza sie do zajetosci.
        """
        bridged = occupied if unknown is None else occupied | unknown
        pad = max_gap + 1                          # zeby zamykanie nie obcinalo koncow
        closed = binary_closing(np.pad(bridged, pad), structure=np.ones(max_gap + 1))[pad:-pad]
        labels, n = label(closed)
        segments = []
        for i in range(1, n + 1):
            idx = np.nonzero(labels == i)[0]
            start, stop = idx[0], idx[-1] + 1
            if stop - start >= min_length and occupied[start:stop].mean() >= min_fill:
                segments.append((start, stop))
        return segments


    @staticmethod
    def pca(image_cut):
        e = np.nan
        t = np.nan
        f = np.nan
        # x = column, y = row; theta is the major-axis angle counter-clockwise
        # from +x, folded into [-pi/2, pi/2) (an orientation has period pi).
        y, x = np.indices(image_cut.shape)
        I = image_cut.clip(min=0)
        It = I.sum()
        if It > 0:
            # second moments about the centroid, not the box centre
            dx = x - np.sum(I * x) / It
            dy = y - np.sum(I * y) / It

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
                t = (np.arctan2(vy, vx) + np.pi / 2) % np.pi - np.pi / 2
        return f,e,t

    @staticmethod
    def adaptive_moments(image_cut, sigma0=2.0, max_iter=30, tol=1e-4):
        """Shape from Gaussian-weighted (adaptive) second moments.

        The elliptical Gaussian window is iterated to match the source
        (Bernstein & Jarvis 2002), so noise and neighbours far from the
        core get little weight and no clipping at zero is needed. For a
        Gaussian PSF the converged moments equal its covariance.

        ``image_cut`` must be background subtracted. Returns
        ``(fwhm, ellipticity, theta)`` like :meth:`pca`, all NaN when the
        iteration does not converge.
        """
        f = e = t = np.nan
        I = np.asarray(image_cut, dtype=float)
        if I.size == 0 or not np.all(np.isfinite(I)):
            return f, e, t
        y, x = np.indices(I.shape)

        cx, cy = FFS.centroid(I)
        if cx is None:
            return f, e, t
        M = np.eye(2) * sigma0 ** 2

        for _ in range(max_iter):
            dx, dy = x - cx, y - cy
            Mi = np.linalg.inv(M)
            r2 = Mi[0, 0] * dx * dx + 2 * Mi[0, 1] * dx * dy + Mi[1, 1] * dy * dy
            wI = np.exp(-0.5 * r2) * I
            S = wI.sum()
            if S <= 0:
                return f, e, t
            cx = cx + np.sum(wI * dx) / S
            cy = cy + np.sum(wI * dy) / S
            # factor 2: for a Gaussian source with window == source the
            # weighted moments are half the true covariance
            M_new = 2.0 * np.array([[np.sum(wI * dx * dx), np.sum(wI * dx * dy)],
                                    [np.sum(wI * dx * dy), np.sum(wI * dy * dy)]]) / S
            if np.linalg.det(M_new) <= 0 or M_new[0, 0] <= 0:
                return f, e, t
            done = np.max(np.abs(M_new - M)) < tol * np.trace(M)
            M = M_new
            if done:
                break
        else:
            return f, e, t

        eigvals, eigvecs = np.linalg.eigh(M)
        b2, a2 = eigvals
        f = 2.355 * (np.sqrt(a2) + np.sqrt(b2)) / 2
        e = 1.0 - np.sqrt(b2 / a2)
        vx, vy = eigvecs[:, 1]
        t = (np.arctan2(vy, vx) + np.pi / 2) % np.pi - np.pi / 2
        return f, e, t

    @staticmethod
    def cpe(image_cut):
        cpe = np.nan
        # central pixel excess
        cx = int(image_cut.shape[0]/2)
        cy = int(image_cut.shape[1]/2)

        if cx > 1 and cy > 1:
            I_std = np.nanstd(image_cut)

            # Zabezpieczenie przed dzieleniem przez zero
            if I_std == 0 or np.isnan(I_std):
                return np.nan

            cx, cy = FFS.centroid(image_cut)
            cx_i = int(round(cx))
            cy_i = int(round(cy))

            I_peak = image_cut[cy_i, cx_i]
            flux = np.sum(image_cut)

            # opcja rekomendowana
            if flux <= 0 or not np.isfinite(flux):
                return np.nan

            cpe = I_peak / flux

        return cpe


#########################################

    @staticmethod
    def fwhm(image_cut):
        if image_cut.size == 0:
            return np.nan, np.nan

        # tylko dodatni sygnał
        data = image_cut.astype(float)
        data[~np.isfinite(data)] = 0
        data = data.clip(min=0)

        if np.sum(data) <= 0:
            return np.nan, np.nan

        # lekkie wygładzenie (kluczowe!)
        data = gaussian_filter(data, sigma=0.5)

        # centroid
        cx, cy = FFS.centroid(data)
        if cx is None or cy is None:
            return np.nan, np.nan

        cx_i = int(round(cx))
        cy_i = int(round(cy))

        # zabezpieczenie zakresu
        ny, nx = data.shape
        if not (0 <= cx_i < nx and 0 <= cy_i < ny):
            return np.nan, np.nan

        # profile
        line_x = data[cy_i, :]
        line_y = data[:, cx_i]

        fwhm_x = FFS.fwhm_1d(line_x)
        fwhm_y = FFS.fwhm_1d(line_y)

        return fwhm_x, fwhm_y

    @staticmethod
    def fwhm_1d(line):
        line = np.asarray(line, dtype=float)

        if not np.any(np.isfinite(line)):
            return np.nan

        # usuń NaNy
        mask = np.isfinite(line)
        if mask.sum() < 3:
            return np.nan
        line = line[mask]



        # znajdź maksimum
        imax = np.argmax(line)

        # peak i half max (bardziej odporne niż max)
        peak = line[imax]
        if peak <= 0:
            return np.nan

        half = 0.5 * peak

        # --- lewa strona ---
        left = imax
        while left > 0 and line[left] > half:
            left -= 1

        if left == 0:
            return np.nan

        # interpolacja
        y1, y2 = line[left], line[left + 1]
        if y2 == y1:
            xl = left
        else:
            xl = left + (half - y1) / (y2 - y1)

        # --- prawa strona ---
        right = imax
        while right < len(line) - 1 and line[right] > half:
            right += 1

        if right == len(line) - 1:
            return np.nan

        y1, y2 = line[right - 1], line[right]
        if y2 == y1:
            xr = right
        else:
            xr = (right - 1) + (half - y1) / (y2 - y1)

        fwhm = xr - xl

        if fwhm <= 0 or not np.isfinite(fwhm):
            return np.nan

        fwhm = np.sqrt(fwhm**2 - (2.355 * 0.5)**2) # poprawka na splot z gaussem sigma = 1

        return fwhm

    ###########################################

################################################

    @staticmethod
    def concentration_index(image_cut, r1=2.0, r2=4.0):
        if image_cut.size == 0:
            return np.nan

        # centroid (wazne!)
        y, x = np.indices(image_cut.shape)
        I = image_cut.clip(min=0)

        flux_total = np.sum(I)
        if flux_total <= 0 or not np.isfinite(flux_total):
            return np.nan

        cy = np.sum(y * I) / flux_total
        cx = np.sum(x * I) / flux_total

        # odleglosci od centroidu
        r = np.sqrt((x - cx) ** 2 + (y - cy) ** 2)

        flux_r1 = np.sum(I[r <= r1])
        flux_r2 = np.sum(I[r <= r2])

        if flux_r2 <= 0:
            return np.nan

        return flux_r1 / flux_r2

    @staticmethod
    def centroid(image_cut):
        I = image_cut.clip(min=0)
        Itot = I.sum()

        if Itot <= 0:
            return None, None

        y, x = np.indices(I.shape)

        cx = np.sum(x * I) / Itot
        cy = np.sum(y * I) / Itot

        return cx, cy

    @staticmethod
    def bin_image(image, bin):
        """Binning bin x bin srednia (ADU w tej samej skali co oryginal; NaN w bloku daje NaN).

        Niepelne bloki na koncu (gdy wymiar nie dzieli sie przez bin) sa odcinane.
        bin=1 zwraca obraz bez zmian (bez kopii).
        """
        if bin == 1:
            return image
        ny, nx = image.shape
        ny, nx = ny // bin * bin, nx // bin * bin
        return image[:ny, :nx].reshape(ny // bin, bin, nx // bin, bin).mean(axis=(1, 3))

    # Przeliczenia obraz zbinowany <-> oryginal (self.bin). Przy bin=1 wszystko jest tozsamoscia.
    # Piksel zbinowany i obejmuje piksele oryginalu i*bin ... i*bin + bin - 1, jego srodek
    # to i*bin + (bin - 1) / 2. Dlugosci (np. FWHM) skaluja sie po prostu przez bin.

    def to_original_coord(self, c):
        """Wspolrzedna (x lub y, skalar lub tablica) z obrazu zbinowanego do oryginalu."""
        return np.asarray(c, dtype=float) * self.bin + (self.bin - 1) / 2.0

    def to_binned_coord(self, c):
        """Wspolrzedna (x lub y, skalar lub tablica) z oryginalu do obrazu zbinowanego."""
        return (np.asarray(c, dtype=float) - (self.bin - 1) / 2.0) / self.bin

    def to_original_mask(self, mask):
        """Maska (lub mapa) z obrazu zbinowanego do wymiarow oryginalu: piksel -> blok bin x bin.

        Wiersze/kolumny odciete przy binowaniu (niepelne bloki) dostaja wartosc sasiedniego bloku.
        """
        if self.bin == 1:
            return mask
        out = np.repeat(np.repeat(mask, self.bin, axis=0), self.bin, axis=1)
        ny, nx = self.image.shape
        return np.pad(out, ((0, ny - out.shape[0]), (0, nx - out.shape[1])), mode="edge")

    def to_binned_mask(self, mask):
        """Maska z oryginalu do obrazu zbinowanego: piksel zbinowany zamaskowany, gdy zamaskowany
        jest ktorykolwiek piksel jego bloku (OR). Niepelne bloki na koncu sa odcinane."""
        if self.bin == 1:
            return mask
        b = self.bin
        ny, nx = mask.shape
        ny, nx = ny // b * b, nx // b * b
        return mask[:ny, :nx].reshape(ny // b, b, nx // b, b).any(axis=(1, 3))

    @staticmethod
    def make_segments(image, n_segments=10, overlap=0):
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

    @staticmethod
    def polysurf(x, y, coeff):
        a0, a1, a2, a3, a4, a5 = coeff
        return a0 + a1 * x + a2 * y + a3 * x ** 2 + a4 * x * y + a5 * y ** 2

    @staticmethod
    def line_filter(image,kernel1_size=3,kernel2_size=7,th1=3,th2=3):

        kernel_l, kernel_r = FFS.line_detection_kernel(kernel1_size)  # tu sie zmienia
        result_l = convolve(image, kernel_l)
        maska_1l = result_l > np.median(result_l) + th1 * mad_std(result_l)
        result_r = convolve(image, kernel_r)
        maska_1r = result_r > np.median(result_r) + th1 * mad_std(result_r)

        kernel_l, kernel_r = FFS.line_detection_kernel(kernel2_size)  # tu sie zmienia
        result_l = convolve(image, kernel_l)
        maska_2l = result_l > np.median(result_l) + th2 * mad_std(result_r) # to nie blad, chcemy miec ta sama skalef
        result_r = convolve(image, kernel_r)
        maska_2r = result_r > np.median(result_r) + th2 * mad_std(result_r)

        maska = (maska_1l & maska_2l) ^ (maska_1r & maska_2r)

        return maska

    @staticmethod
    def line_detection_kernel(kernel_half_size=5):
        size = 2 * kernel_half_size + 1
        center = kernel_half_size
        y, x = np.ogrid[:size, :size]
        distance = np.sqrt((x - center) ** 2 + (y - center) ** 2)
        inner = kernel_half_size - 1 / 2
        outer = kernel_half_size + 1 / 2

        # The four ring pixels on the axes are split between the halves: the
        # horizontal-axis pair (y == center) goes to "left", the vertical-axis
        # pair (x == center) to "right", so each kernel sums to zero and no
        # line direction gives a zero response.
        left = ((x < center) & (y <= center)) | ((x > center) & (y >= center))
        right = ((x <= center) & (y > center)) | ((x >= center) & (y < center))

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

    @staticmethod
    def laplace_kernel():
        kernel = [[0,-1,0],[-1,4,-1],[0,-1,0]]
        return kernel

    @staticmethod
    def gauss2d_kernel(size, sigma):
        kernel = np.fromfunction(lambda x, y: (1 / (2 * np.pi * sigma ** 2)) * np.exp(
            -((x - (size - 1) / 2) ** 2 + (y - (size - 1) / 2) ** 2) / (2 * sigma ** 2)), (size, size))
        return kernel / np.sum(kernel)

    @staticmethod
    def _extrapolate(coords, values, edge):
        """Value at edge linearly extrapolated from the two finite values nearest to it
        (one finite value - that value; none - NaN). coords sorted ascending."""
        ok = np.isfinite(values)
        c, v = coords[ok], values[ok]
        if len(c) == 0:
            return np.nan
        if len(c) == 1:
            return v[0]
        if edge <= c[0]:
            (c1, v1), (c2, v2) = (c[0], v[0]), (c[1], v[1])
        else:
            (c1, v1), (c2, v2) = (c[-1], v[-1]), (c[-2], v[-2])
        return v1 + (edge - c1) * (v2 - v1) / (c2 - c1)

    @staticmethod
    def _column_deviation(resid, block, window):
        # mediana kazdej kolumny w blokach po block wierszy minus mediana window sasiednich kolumn
        ny, nx = resid.shape
        nb = max(1, ny // block)
        rows = np.arange(ny) * nb // ny                     # numer bloku kazdego wiersza
        prof = np.array([np.nanmedian(resid[rows == i], axis=0) for i in range(nb)])
        dev = prof - median_filter(prof, size=(1, window), mode="nearest")
        return dev[rows]

    @staticmethod
    def _fits_text(text):
        """Text for a FITS header: typographic characters -> ASCII, the rest of non-ASCII dropped."""
        text = text.translate(str.maketrans({"–": "-", "−": "-", "×": "x"}))
        return text.encode("ascii", "ignore").decode()

    @staticmethod
    def line_to_mask(shape, rho, theta, half_width=0.5, x0=None, y0=None, x1=None, y1=None):
        """Maska pikseli linii x*cos(theta) + y*sin(theta) = rho (konwencja find_lines, x = kolumna, y = wiersz).

        Zaznacza piksele w odleglosci <= half_width od linii. Gdy podane sa konce odcinka
        (x0, y0, x1, y1), tylko odcinek (z zapasem half_width na koncach), bez nich cala linia
        przez kadr. Zwraca maske bool o ksztalcie shape; np.nonzero(maska) daje (y, x) pikseli.
        """
        mask = np.zeros(shape, dtype=bool)
        ny, nx = shape
        c, s = np.cos(theta), np.sin(theta)
        if x0 is None:
            ends = FFS._frame_crossings(shape, rho, theta)
            if ends is None:
                return mask                    # linia nie przechodzi przez kadr
            x0, y0, x1, y1 = ends
        pad = half_width + 1
        xa, xb = int(max(0, min(x0, x1) - pad)), int(min(nx, max(x0, x1) + pad + 1))
        ya, yb = int(max(0, min(y0, y1) - pad)), int(min(ny, max(y0, y1) + pad + 1))
        yy, xx = np.mgrid[ya:yb, xa:xb]
        dist = np.abs(xx * c + yy * s - rho)            # odleglosc od linii
        pos = -xx * s + yy * c                         # polozenie wzdluz linii
        p0, p1 = sorted([-x0 * s + y0 * c, -x1 * s + y1 * c])
        mask[ya:yb, xa:xb] = (dist <= half_width) & (pos >= p0 - half_width) & (pos <= p1 + half_width)
        return mask

    @staticmethod
    def _frame_crossings(shape, rho, theta):
        # konce odcinka, w ktorym linia przecina kadr [0, nx-1] x [0, ny-1]; None gdy go omija
        ny, nx = shape
        c, s = np.cos(theta), np.sin(theta)
        pts = []
        if abs(s) > 1e-12:                     # przeciecia z lewa i prawa krawedzia
            for x in (0, nx - 1):
                y = (rho - x * c) / s
                if 0 <= y <= ny - 1:
                    pts.append((x, y))
        if abs(c) > 1e-12:                     # przeciecia z gorna i dolna krawedzia
            for y in (0, ny - 1):
                x = (rho - y * s) / c
                if 0 <= x <= nx - 1:
                    pts.append((x, y))
        if not pts:
            return None
        # najdalsza para punktow (w rogach te same punkty pojawiaja sie dwa razy)
        (xa, ya), (xb, yb) = max(((p, q) for p in pts for q in pts),
                                 key=lambda pq: (pq[0][0] - pq[1][0]) ** 2 + (pq[0][1] - pq[1][1]) ** 2)
        return xa, ya, xb, yb
