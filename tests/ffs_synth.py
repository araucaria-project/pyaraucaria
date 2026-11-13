"""Synthetic CCD frames with known truth, for testing ``pyaraucaria.ffs``.

Conventions (the ones FFS is expected to follow):

- ``image[y, x]`` -- first index is the row (y), second the column (x), 0-based.
- A star at ``(x, y)`` has its centre at column ``x``, row ``y``.
- ``theta`` is the position angle of the PSF major axis, in radians,
  counter-clockwise from the +x axis (towards +y).
- ``fwhm`` is the major-axis FWHM; minor axis is ``fwhm * (1 - ellipticity)``.

Every artifact injected into the frame is also recorded as a boolean mask in
``SyntheticFrame.truth`` so that future mask-producing code can be scored
against it.
"""

from dataclasses import dataclass, field

import numpy as np

FWHM_TO_SIGMA = 1.0 / (2.0 * np.sqrt(2.0 * np.log(2.0)))


@dataclass
class SyntheticFrame:
    image: np.ndarray
    stars: list
    truth: dict = field(default_factory=dict)
    params: dict = field(default_factory=dict)


def gaussian_star(shape, x, y, flux, fwhm, ellipticity=0.0, theta=0.0, radius=None):
    """Return ``(slices, stamp)`` of an elliptical Gaussian with total ``flux``."""
    sig_a = fwhm * FWHM_TO_SIGMA
    sig_b = sig_a * (1.0 - ellipticity)
    if radius is None:
        radius = int(np.ceil(5 * sig_a)) + 1
    ny, nx = shape
    x0, x1 = max(0, int(x) - radius), min(nx, int(x) + radius + 1)
    y0, y1 = max(0, int(y) - radius), min(ny, int(y) + radius + 1)
    yy, xx = np.mgrid[y0:y1, x0:x1]
    dx = xx - x
    dy = yy - y
    c, s = np.cos(theta), np.sin(theta)
    u = dx * c + dy * s       # along major axis
    v = -dx * s + dy * c      # along minor axis
    stamp = np.exp(-0.5 * ((u / sig_a) ** 2 + (v / sig_b) ** 2))
    stamp *= flux / (2.0 * np.pi * sig_a * sig_b)
    return (slice(y0, y1), slice(x0, x1)), stamp


def _segment_distance(xx, yy, x0, y0, x1, y1):
    px, py = x1 - x0, y1 - y0
    norm = px * px + py * py
    t = ((xx - x0) * px + (yy - y0) * py) / norm
    t = np.clip(t, 0.0, 1.0)
    return np.hypot(xx - (x0 + t * px), yy - (y0 + t * py))


def make_frame(shape=(512, 512), sky=1000.0, sky_gradient=(0.0, 0.0),
               sky_curvature=0.0, gain=1.0, read_noise=5.0, stars=(),
               hot_pixels=(), bad_columns=(), lines=(), saturation=None,
               bleed=True, dtype=np.float64, seed=0):
    """Build a synthetic frame.

    Parameters
    ----------
    shape : (ny, nx)
    sky : float
        Sky level in ADU at the frame centre.
    sky_gradient : (gx, gy)
        Linear sky gradient in ADU per pixel along x and y.
    sky_curvature : float
        Quadratic term ``c * r^2`` (r from the centre, ADU / px^2).
    gain : float
        e-/ADU. Poisson noise is applied in electrons.
    read_noise : float
        Read noise in electrons.
    stars : iterable of dict
        Keys ``x, y, flux`` (ADU, total), optional ``fwhm`` (default 4),
        ``ellipticity`` (0), ``theta`` (0).
    hot_pixels : iterable of (x, y, adu)
        Added on top of everything, after noise.
    bad_columns : iterable of dict
        Keys ``x``, ``kind`` ("hot" adds ``value`` ADU, "dead" multiplies the
        column by ``value``), optional ``y0``/``y1`` row range.
    lines : iterable of dict
        Satellite trails: ``x0, y0, x1, y1`` (segment ends), ``flux`` (peak ADU
        per pixel), ``width`` (FWHM across the trail, px).
    saturation : float or None
        Clip level in ADU. With ``bleed=True`` excess charge of each column is
        spread up and down from the saturated run.
    dtype : numpy dtype
        Output dtype; integer dtypes are rounded and clipped to range.
    seed : int
        Seed for the noise generator; same seed -> identical frame.
    """
    rng = np.random.default_rng(seed)
    ny, nx = shape
    yy, xx = np.mgrid[0:ny, 0:nx].astype(float)
    xc, yc = (nx - 1) / 2.0, (ny - 1) / 2.0

    model = (sky + sky_gradient[0] * (xx - xc) + sky_gradient[1] * (yy - yc)
             + sky_curvature * ((xx - xc) ** 2 + (yy - yc) ** 2))
    sky_model = model.copy()

    star_list = []
    for st in stars:
        st = {"fwhm": 4.0, "ellipticity": 0.0, "theta": 0.0, **st}
        sl, stamp = gaussian_star(shape, st["x"], st["y"], st["flux"], st["fwhm"],
                                  st["ellipticity"], st["theta"])
        model[sl] += stamp
        st["peak"] = float(stamp.max())
        star_list.append(st)

    mask_line = np.zeros(shape, dtype=bool)
    for ln in lines:
        sigma = ln.get("width", 2.0) * FWHM_TO_SIGMA
        d = _segment_distance(xx, yy, ln["x0"], ln["y0"], ln["x1"], ln["y1"])
        model += ln["flux"] * np.exp(-0.5 * (d / sigma) ** 2)
        mask_line |= d <= 3 * sigma

    electrons = np.clip(model * gain, 0, None)
    data = rng.poisson(electrons).astype(float)
    data += rng.normal(0.0, read_noise, size=shape)
    data /= gain

    mask_bad_column = np.zeros(shape, dtype=bool)
    for bc in bad_columns:
        y0, y1 = bc.get("y0", 0), bc.get("y1", ny)
        if bc.get("kind", "hot") == "dead":
            data[y0:y1, bc["x"]] *= bc.get("value", 0.0)
        else:
            data[y0:y1, bc["x"]] += bc.get("value", 500.0)
        mask_bad_column[y0:y1, bc["x"]] = True

    mask_hot = np.zeros(shape, dtype=bool)
    for x, y, adu in hot_pixels:
        data[y, x] += adu
        mask_hot[y, x] = True

    mask_saturated = np.zeros(shape, dtype=bool)
    if saturation is not None:
        over = data > saturation
        if bleed and over.any():
            for x in np.nonzero(over.any(axis=0))[0]:
                col = data[:, x]
                rows = np.nonzero(col > saturation)[0]
                excess = float(np.sum(col[rows] - saturation))
                extra = int(excess // saturation)
                lo = max(0, rows.min() - extra // 2)
                hi = min(ny, rows.max() + 1 + extra - extra // 2)
                col[lo:hi] = saturation
        mask_saturated = data >= saturation
        data = np.minimum(data, saturation)

    if np.issubdtype(np.dtype(dtype), np.integer):
        info = np.iinfo(dtype)
        data = np.clip(np.round(data), info.min, info.max)
    data = data.astype(dtype)

    truth = {
        "sky_model": sky_model,
        "mask_hot": mask_hot,
        "mask_bad_column": mask_bad_column,
        "mask_line": mask_line,
        "mask_saturated": mask_saturated,
        "lines": [dict(ln) for ln in lines],
    }
    params = dict(shape=shape, sky=sky, sky_gradient=sky_gradient,
                  sky_curvature=sky_curvature, gain=gain, read_noise=read_noise,
                  saturation=saturation, seed=seed)
    return SyntheticFrame(image=data, stars=star_list, truth=truth, params=params)


def star_field(n=40, shape=(512, 512), flux_range=(2e3, 2e5), fwhm=4.0,
               ellipticity=0.0, theta=0.0, margin=20, min_sep=15, seed=1):
    """Random, well separated stars with log-uniform flux."""
    rng = np.random.default_rng(seed)
    ny, nx = shape
    out = []
    tries = 0
    while len(out) < n and tries < 100 * n:
        tries += 1
        x = rng.uniform(margin, nx - margin)
        y = rng.uniform(margin, ny - margin)
        if any(np.hypot(x - s["x"], y - s["y"]) < min_sep for s in out):
            continue
        flux = float(np.exp(rng.uniform(*np.log(flux_range))))
        out.append(dict(x=x, y=y, flux=flux, fwhm=fwhm,
                        ellipticity=ellipticity, theta=theta))
    return out


# --- Named scenarios shared by regression and truth tests -------------------

def scenario(name):
    """Return a deterministic ``SyntheticFrame`` by name."""
    if name == "typical":
        stars = star_field(n=40, seed=1)
        # Two bright stars that will saturate and bleed.
        stars += [dict(x=150.3, y=380.7, flux=3e6, fwhm=4.0),
                  dict(x=400.6, y=120.2, flux=2e6, fwhm=4.0)]
        return make_frame(
            stars=stars, sky=1000.0, sky_gradient=(0.3, -0.2),
            hot_pixels=[(37, 211, 8000), (300, 45, 20000), (455, 470, 3000)],
            bad_columns=[dict(x=260, kind="hot", value=400.0)],
            lines=[dict(x0=0, y0=60, x1=511, y1=300, flux=300.0, width=2.5)],
            saturation=45000, seed=2)
    if name == "typical_uint16":
        f = scenario("typical")
        f.image = np.clip(np.round(f.image), 0, 65535).astype(np.uint16)
        return f
    if name == "clean":
        return make_frame(stars=star_field(n=30, flux_range=(2e4, 2e5), seed=3),
                          sky=500.0, seed=4)
    if name == "elliptical":
        return make_frame(stars=star_field(n=30, flux_range=(2e4, 2e5), fwhm=5.0,
                                           ellipticity=0.3, theta=np.deg2rad(30),
                                           seed=5),
                          sky=500.0, seed=6)
    if name == "sparse":
        return make_frame(stars=star_field(n=3, flux_range=(5e4, 1e5), seed=7),
                          sky=500.0, seed=8)
    if name == "empty":
        return make_frame(sky=500.0, seed=9)
    if name == "gradient":
        return make_frame(shape=(400, 600), sky=1000.0, sky_gradient=(0.5, 0.25),
                          sky_curvature=1e-3, seed=10)
    raise KeyError(name)


SCENARIOS = ["typical", "typical_uint16", "clean", "elliptical", "sparse",
             "empty", "gradient"]
