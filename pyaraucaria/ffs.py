import numpy as np
from scipy.ndimage import maximum_filter, median_filter
from scipy.ndimage import convolve
from scipy.ndimage import gaussian_filter
from scipy.ndimage import binary_closing, label
from scipy.interpolate import LinearNDInterpolator, NearestNDInterpolator

from astropy.stats import mad_std
from astropy.table import Table



class FFS:

    def __init__(self, image, gain=1., rn_noise=0., saturation = 50000):
        self.image = image
        self.gain = float(gain)
        self.rn_noise = float(rn_noise)
        self.saturation = saturation

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
        self.masks = {}                # maski logiczne w rozmiarze obrazu, osobno dla kazdej metody
        self.masks_description = {}
        self.exclude = set()           # nazwy masek wykluczajacych (zle piksele), patrz bad_mask()

        nonfinite = ~np.isfinite(self.image)
        if nonfinite.any():            # NaN/inf na wejsciu (np. po redukcji) traktujemy jak zle piksele
            self.masks["nonfinite"] = nonfinite
            self.masks_description["nonfinite"] = "Non-finite input pixels (NaN or inf)"
            self.exclude.add("nonfinite")

        self.maps = {}
        self.maps_description = {}

    def bad_mask(self, ignore=()):
        """Suma (OR) masek wykluczajacych, ktorych nazwy sa w self.exclude (poza nazwami z ignore)."""
        bad = np.zeros(self.image.shape, dtype=bool)
        for name in self.exclude:
            if name not in ignore:
                bad |= self.masks[name]
        return bad

    def _background(self):
        """Poziom tla i sigma: (maps["sky"], sigma residuum) jesli jest mapa tla, inaczej (median, q_sigma).

        Sigma z residuum image - sky liczona z kwantyli dobrych pikseli: q_sigma calej klatki
        zawiera tez rozrzut samego gradientu tla. Wymaga mk_stats().
        """
        if "sky" not in self.maps:
            return self.median, self.q_sigma, "median"
        bad = self.bad_mask()
        residual = self.image - self.maps["sky"]
        good = residual[~bad] if bad.sum() < bad.size else residual.ravel()
        q16, q84 = np.percentile(good, [15.9, 84.1])
        return self.maps["sky"], (q84 - q16) / 2.0, "sky map"

    def mk_stats(self):
        # statystyki z dobrych pikseli (poza bad_mask); min/max ze wszystkich, zeby bylo widac saturacje
        bad = self.bad_mask()
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

        self.mk_threshold_mask()

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

    def mk_threshold_mask(self, nsigma=3.0):
        """Maska pikseli powyzej tla: image > tlo + nsigma * sigma.

        Tlo to self.maps["sky"], jesli juz policzone (sky_map); wtedy sigma liczona jest
        z kwantyli residuum image - sky. W przeciwnym razie mediana i q_sigma z mk_stats().
        """
        level, sigma, level_name = self._background()
        self.masks["threshold"] = (self.image > level + nsigma * sigma) & ~self.bad_mask()
        self.masks_description["threshold"] = (
            f"Pixels above {level_name} + {nsigma} * sigma (sigma = {sigma:.2f} ADU, quantile-based; "
            f"excluding bad_mask; input mask for line detection)"
        )

    def mk_saturation_mask(self, saturation=None):
        """Maska pikseli nasyconych: image >= saturation (domyslnie self.saturation)."""
        if saturation is None:
            saturation = self.saturation
        self.masks["saturation"] = self.image >= saturation
        self.exclude.add("saturation")
        self.masks_description["saturation"] = f"Saturated pixels (value >= {saturation} ADU)"

    def mk_columns_map(self, block=64, window=15):
        """Mapy odchylen kolumn i wierszy (klatki kalibracyjne: zera, darki, flaty, bez gwiazd).

        Od obrazu odejmowana jest mapa tla (sky_map, liczona jesli jej nie ma; liczbe segmentow
        mozna ustawic wczesniej przez sky_map(n_segments=...)). maps["columns"] to roznica mediany
        kolumny w bloku block wierszy i mediany window sasiednich kolumn, w ADU, w rozmiarze obrazu.
        maps["rows"] tak samo dla wierszy. Maski z progiem: mk_bad_columns_mask().
        """
        if "sky" not in self.maps:
            self.sky_map()
        resid = np.where(self.bad_mask(), np.nan, self.image - self.maps["sky"])
        self.maps["columns"] = FFS._column_deviation(resid, block, window)
        self.maps["rows"] = FFS._column_deviation(resid.T, block, window).T
        for name, what in (("columns", "column"), ("rows", "row")):
            self.maps_description[name] = (
                f"{what.capitalize()} median in blocks of {block} px minus the median of {window} neighbouring "
                f"{what}s, sky map subtracted (ADU)"
            )

    def mk_bad_columns_mask(self, threshold):
        """Maski zlych kolumn i wierszy: |maps["columns"]| > threshold, |maps["rows"]| > threshold (ADU).

        Wymaga mk_columns_map(). Obie maski sa dopisywane do self.exclude.
        """
        for name, what in (("columns", "column"), ("rows", "row")):
            self.masks["bad_" + name] = np.abs(self.maps[name]) > threshold
            self.masks_description["bad_" + name] = (
                f"Bad {what}s: |maps['{name}']| > {threshold} ADU (mk_bad_columns_mask)"
            )
            self.exclude.add("bad_" + name)

    def mk_pixels_map(self, size=5):
        """Mapa odchylen pojedynczych pikseli: image - median_filter(image, size x size), w ADU.

        Dla klatek kalibracyjnych (najlepiej master; na pojedynczym darku hot piksel i CR wygladaja
        tak samo). Maski z progiem: mk_hot_pixels_mask(), mk_cold_pixels_mask().
        """
        image = self.image.astype(float)
        self.maps["pixels"] = image - median_filter(image, size=size)
        self.maps_description["pixels"] = f"Pixel value minus the median of its {size}x{size} neighbourhood (ADU)"

    def mk_hot_pixels_mask(self, threshold):
        """Maska hot pikseli: maps["pixels"] > threshold (ADU). Wymaga mk_pixels_map()."""
        self.masks["hot_pixels"] = self.maps["pixels"] > threshold
        self.masks_description["hot_pixels"] = f"Hot pixels: maps['pixels'] > {threshold} ADU (mk_hot_pixels_mask)"
        self.exclude.add("hot_pixels")

    def mk_cold_pixels_mask(self, threshold):
        """Maska cold pikseli: maps["pixels"] < -threshold (ADU). Wymaga mk_pixels_map()."""
        self.masks["cold_pixels"] = self.maps["pixels"] < -threshold
        self.masks_description["cold_pixels"] = f"Cold pixels: maps['pixels'] < -{threshold} ADU (mk_cold_pixels_mask)"
        self.exclude.add("cold_pixels")

    def find_stars(self, threshold=5, method="sigma quantile", fwhm=10,
                   min_smoothed_sigma=None, rank_by="raw",
                   max_concentration=None, aperture_radius=None,
                   min_pixels_above_threshold=None):
        """Detect star peaks in ``self.image``.

        Parameters
        ----------
        threshold : float
            Detection threshold in units of background sigma (raw image).
        method : str
            How background sigma is estimated — see class docs.
        fwhm : float
            Adopted PSF FWHM (pixels). Drives the Gaussian kernel sigma
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
            Half-width in pixels of the square aperture used for both
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

        # poziom tla: mapa tla jesli policzona (sky_map), inaczej mediana
        level, q_sigma, _ = self._background()

        if self.fs_method == "rms Poisson":
            self.fs_sigma = self.noise
        elif self.fs_method == "rms":
            self.fs_sigma = self.rms
        elif self.fs_method == "sigma quantile":
            self.fs_sigma = q_sigma
        else:
            raise ValueError(f"Invalid method type {self.fs_method}")

        # Cast to float first: gaussian_filter preserves the input dtype, so on
        # integer FITS data (uint16/int16) the smoothed image — and every test
        # built on it (min_smoothed_sigma, rank_by="smoothed") — would be
        # truncated to integers.
        image = self.image.astype(float)
        # zle piksele (bad_mask, NaN) wypelniamy tlem, zeby nie byly kandydatami i nie psuly
        # wygladzania; saturacji nie, bo nasycone gwiazdy maja zostac w coo/adu
        bad = self.bad_mask(ignore=("saturation",))
        image[bad] = level[bad] if np.ndim(level) else level

        mask1 = image > level + self.fs_threshold * self.fs_sigma
        data2 = gaussian_filter(image, sigma=self.fs_kernel_sigma)
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
        val = self.image[mask]

        # Per-candidate aperture statistics (background-subtracted) —
        # needed for max_concentration filter and rank_by="aperture".
        # Computed once for both, since both are cheap (~50 ops/cand).
        ap_radius = self.fs_aperture_radius
        H, W = self.image.shape
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
            min_px = int(min_pixels_above_threshold)
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
        self.mk_saturation_mask()
        self.mk_stats()
        self.sky_map()
        self.mk_threshold_mask()
        self.find_lines()
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

        mk = np.ones_like(self.fwhm, dtype=bool)

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

    def star_info(self, box=10, N_stars=None):

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

        if len(self.coo) < 5:
            return

        # Measurements of accepted stars only; `keep` holds their indices in
        # self.coo so every output row stays aligned with its source.
        bad = self.bad_mask()
        keep = []
        rows = []

        for i, (y, x) in enumerate(self.coo):

            if len(keep) >= N_stars:
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

            keep.append(i)
            rows.append((box_mag, bkg, e, t, star_fwhm, fx, fy, cpe, shape, ci))

        cols = np.array(rows, dtype=float).reshape(-1, 10).T
        (self.box_mag, self.bkg, self.ellipticity, self.theta, self.fwhm,
         self.fwhm_x, self.fwhm_y, self.cpe, self.shape, self.ci) = cols
        keep = np.array(keep, dtype=int)
        self.coo = self.coo[keep]
        self.adu = self.adu[keep]

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

    def sky_map(self,n_segments=10):
        self.sky_gradient(n_segments=n_segments)

        points = np.column_stack((self.stats["sky"]["surface_x"], self.stats["sky"]["surface_y"]))
        values = self.stats["sky"]["surface_val"]

        ny, nx = self.image.shape
        yy, xx = np.mgrid[0:ny, 0:nx]

        sky_map = LinearNDInterpolator(points, values)(xx, yy)

        # poza otoczka wypukla srodkow segmentow (brzegi kadru) LinearNDInterpolator daje NaN
        outside = np.isnan(sky_map)
        if outside.any():
            sky_map[outside] = NearestNDInterpolator(points, values)(xx[outside], yy[outside])

        self.maps["sky"] = sky_map
        self.maps_description["sky"] = (
            "Sky background map: linear interpolation (Delaunay triangulation) of segment medians; "
            "outside the convex hull of segment centres filled with the nearest segment median"
        )

    def sky_gradient(self,n_segments=10, min_good_fraction=0.5):
        image = self.image
        segments = FFS.make_segments(image, n_segments=n_segments)
        bad_segments = FFS.make_segments(self.bad_mask(), n_segments=n_segments)

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
            bk.append(FFS.polysurf(xi, yi, coeff))
            le.append(FFS.polysurf(x0, yi, coeff))
            re.append(FFS.polysurf(xk, yi, coeff))
            ue.append(FFS.polysurf(xi, yk, coeff))
            de.append(FFS.polysurf(xi, y0, coeff))

        max_amplitude = max(bk) - min(bk)
        # co to po co to
        # max(le) - min(re)
        frame_gradient = max(max((max(le) - min(re)), (max(re) - min(le))),
                             max((max(ue) - min(de)), (max(de) - min(ue))))

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
                "Maximum amplitude of the background signal across the image (in ADU)"
            ),

            "bkg_frame_gradient": (
                "Background gradient across the image frame (in ADU)"
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
                   max_candidates=100, steps=None, max_shift=1.0):
        """Detekcja sladow (linii) na podstawie transformaty Hougha maski self.masks["threshold"].

        Kandydaci to komorki akumulatora z najwieksza liczba glosow. Kandydat jest
        linia, jesli wzdluz niej piksele maski tworza ciagly odcinek: przerwy
        <= max_gap px, dlugosc >= 10 * fwhm, zajetosc >= min_fill. Po przyjeciu
        linii jej piksele (w odleglosci <= half_width) i ich glosy sa usuwane.
        steps i max_shift sa przekazywane do hough_transform().
        """
        self.hough_transform(steps=steps, max_shift=max_shift)

        rh, theta = self.rh, self.hough_theta
        rh0 = -int(rh[0])                     # koszyk k <-> rho = k - rh0
        cos_t, sin_t = np.cos(theta), np.sin(theta)
        ys, xs = np.nonzero(self.masks["threshold"])
        xs = xs.astype(float)
        ys = ys.astype(float)

        # znajdujemy maksimum, jak jest powyzej th, a nastepnie usuwamy z glosowania piksele ktore sa w odleglosci
        # half_width od lini

        acc = self.accumulator.copy()
        work = self.masks["threshold"].copy()       # maska, z ktorej usuwamy piksele przyjetych linii
        removed = self.bad_mask()      # piksele "nieznane": zle (bad_mask) i usuniete razem z przyjetymi liniami
        lines_rho, lines_theta, lines_val = [], [], []
        lines_x0, lines_y0, lines_x1, lines_y1 = [], [], [], []

        for _ in range(max_candidates):
            ri_max, ti_max = np.unravel_index(np.argmax(acc), acc.shape)
            val = acc[ri_max, ti_max]
            if val < min_length:
                break
            rho, t = rh[ri_max], theta[ti_max]

            # czy wzdluz linii piksele maski tworza ciagly odcinek?
            pos, occupied = FFS._line_profile(work, rho, t)
            _, unknown = FFS._line_profile(removed, rho, t)   # skrzyzowania z juz znalezionymi liniami
            segments = FFS._line_segments(occupied, max_gap, 10 * fwhm, min_fill, unknown)
            if not segments:
                # odrzucony: wylaczamy te komorke i jej sasiadow, piksele zostaja
                acc[max(0, ri_max - 2):ri_max + 3, max(0, ti_max - 1):ti_max + 2] = -1
                continue

            # przyjety: zapisujemy linie i konce najdluzszego odcinka
            start, stop = max(segments, key=lambda seg: seg[1] - seg[0])
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
            on_line = np.abs(xs * c + ys * s - rho) <= half_width
            FFS._hough_votes(acc, xs[on_line], ys[on_line], cos_t, sin_t, rh0, weight=-1)
            work[ys[on_line].astype(int), xs[on_line].astype(int)] = False
            removed[ys[on_line].astype(int), xs[on_line].astype(int)] = True
            xs, ys = xs[~on_line], ys[~on_line]

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
        self.exclude.add("lines")

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

    def hough_transform(self, steps=None, max_shift=1.0):
        """Transformata Hougha maski self.masks["threshold"].

        Zapisuje self.accumulator (koszyki rho x katy theta), self.rh (rho koszykow,
        krok 1 px) i self.hough_theta. steps=None dobiera liczbe katow tak, zeby
        najdluzsza linia (przekatna) nie odjechala na koncach o wiecej niz max_shift px.
        """
        # Ensure that the mask has been computed or set before use
        if "threshold" not in self.masks:
            raise RuntimeError(
                "self.masks['threshold'] is not set. Call 'mk_stats()' or "
                "'mk_threshold_mask()' first, or set it to a boolean mask array before "
                "calling 'hough_transform()'."
            )
        ys, xs = np.nonzero(self.masks["threshold"])  # piksele ktore biora udzial w zabawie
        xs = xs.astype(float)
        ys = ys.astype(float)

        ny, nx = self.image.shape
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
    def _column_deviation(resid, block, window):
        # mediana kazdej kolumny w blokach po block wierszy minus mediana window sasiednich kolumn
        ny, nx = resid.shape
        nb = max(1, ny // block)
        rows = np.arange(ny) * nb // ny                     # numer bloku kazdego wiersza
        prof = np.array([np.nanmedian(resid[rows == i], axis=0) for i in range(nb)])
        dev = prof - median_filter(prof, size=(1, window), mode="nearest")
        return dev[rows]

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
