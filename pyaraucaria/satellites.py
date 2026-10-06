"""Detection of satellite trails and other straight linear features.

Pipeline (detection runs on a binned, background-flattened copy):

1. Bin the image, subtract a smooth local background (block medians,
   bilinear), and threshold a lightly smoothed residual at ``k_input`` sigma.
   Stars are *not* removed -- a bright trail is itself detected as a chain of
   "stars". Each pixel votes with weight 1 and the rho bin is 1 px wide, so a
   star contributes only its chord to any line.
2. Hough transform of the thresholded pixels. A cell is a candidate when its
   votes exceed the noise expectation ``p * L`` (``p`` -- fraction of pixels
   above threshold, ``L`` -- chord length of that line through the frame) by
   ``z_min`` binomial sigmas.
3. Each candidate is refined locally and verified on the residual image:
   the band along the line is cut into segments, each giving the band mean
   minus the local background from flanking bands. Noise per segment is
   measured empirically on parallel offset lines. The extent is the
   maximum-sum run of capped segment scores (a single star cannot bridge a
   long gap), and the trail is accepted on the *median* segment signal --
   insensitive to stars that happen to lie on the line.
4. Accepted lines are refined at full resolution (centre, width, peak), their
   votes and pixels are removed and the search continues, so one trail is
   reported once. Trails reaching the excluded border are extended to the
   frame edge.

Lines aligned with detector columns or rows within ``axis_tol`` are reported
with ``kind`` "column"/"row" (bad columns, bleed trails, bad rows) and are not
included in the satellite mask.

Conventions: ``image[y, x]``; a line is ``x cos(theta) + y sin(theta) = rho``
with theta in [-pi/2, pi/2), i.e. theta is the direction of the line normal,
counter-clockwise from +x. ``t = -x sin(theta) + y cos(theta)`` is the
position along the line.
"""

import numpy as np
from astropy.stats import mad_std
from astropy.table import Table
from scipy.ndimage import binary_dilation, find_objects, gaussian_filter, label

COLUMNS = ["x0", "y0", "x1", "y1", "rho", "theta", "length", "width", "peak",
           "peak_snr", "mask_halfwidth", "significance", "fraction", "kind"]

DESCRIPTIONS = {
    "x0": "X of the first end of the trail segment (pixels, 0-based)",
    "y0": "Y of the first end of the trail segment (pixels, 0-based)",
    "x1": "X of the second end of the trail segment (pixels, 0-based)",
    "y1": "Y of the second end of the trail segment (pixels, 0-based)",
    "rho": "Distance of the line from the origin: x cos(theta) + y sin(theta) = rho (pixels)",
    "theta": "Angle of the line normal, counter-clockwise from +x, in [-pi/2, pi/2) (radians)",
    "length": "Length of the detected segment (pixels)",
    "width": "FWHM of the trail across its direction (pixels)",
    "peak": "Peak brightness of the trail above local background (ADU per pixel)",
    "peak_snr": "Peak brightness in units of the local pixel noise",
    "mask_halfwidth": "Half-width of the masked band around the trail (pixels)",
    "significance": "Significance of the median segment signal (empirical sigma)",
    "fraction": "Fraction of segments within the extent with significant signal",
    "kind": ("'satellite'; 'column'/'row' for features aligned with the detector axes; "
             "'pattern' for several faint parallel bands (flat/detector pattern). "
             "Only 'satellite' rows are masked"),
}


# --- small helpers ---------------------------------------------------------

def _auto_bin(shape):
    n = max(shape)
    return 1 if n <= 1024 else 2 if n <= 2048 else 4


def _bin(image, b):
    if b == 1:
        return image.astype(float)
    ny, nx = image.shape
    ny2, nx2 = ny - ny % b, nx - nx % b
    return image[:ny2, :nx2].astype(float).reshape(ny2 // b, b, nx2 // b, b).mean(axis=(1, 3))


def _lerp_axis(centres, values, coords, axis):
    """Piecewise-linear interpolation along ``axis`` with linear extrapolation."""
    i = np.clip(np.searchsorted(centres, coords) - 1, 0, len(centres) - 2)
    w = (coords - centres[i]) / (centres[i + 1] - centres[i])
    lo = np.take(values, i, axis=axis)
    hi = np.take(values, i + 1, axis=axis)
    shape = [1] * values.ndim
    shape[axis] = -1
    w = w.reshape(shape)
    return lo * (1 - w) + hi * w


def _background(img, box):
    """Smooth background: block medians on a coarse grid, bilinear.

    Outside the outermost block centres the surface is extrapolated linearly
    -- a flat extension would leave a gradient-dependent offset along the
    frame edges, which Hough then takes for lines parallel to the edges.
    """
    ny, nx = img.shape
    gy, gx = max(1, ny // box), max(1, nx // box)
    if gy < 2 or gx < 2:
        return np.full(img.shape, np.median(img))
    by, bx = ny // gy, nx // gx
    blocks = img[:gy * by, :gx * bx].reshape(gy, by, gx, bx)
    grid = np.median(blocks.transpose(0, 2, 1, 3).reshape(gy, gx, -1), axis=2)
    yc = (np.arange(gy) + 0.5) * by - 0.5
    xc = (np.arange(gx) + 0.5) * bx - 0.5
    rows = _lerp_axis(xc, grid, np.arange(nx, dtype=float), axis=1)     # (gy, nx)
    return _lerp_axis(yc, rows, np.arange(ny, dtype=float), axis=0)      # (ny, nx)


def _hough(ys, xs, ws, cos_t, sin_t, rho_max, chunk=4000):
    """Weighted vote accumulator, shape (2*rho_max+1, n_theta), rho bins of 1 px."""
    T = len(cos_t)
    R = 2 * rho_max + 1
    acc = np.zeros(R * T)
    cols = np.arange(T)
    for i in range(0, len(xs), chunk):
        r = np.rint(xs[i:i + chunk, None] * cos_t + ys[i:i + chunk, None] * sin_t).astype(np.int64)
        w = np.broadcast_to(ws[i:i + chunk, None], r.shape)
        acc += np.bincount(((r + rho_max) * T + cols).ravel(), weights=w.ravel(), minlength=R * T)
    return acc.reshape(R, T)


def _component_weights(hot, area_cap, max_fill=0.25):
    """Per-pixel vote weights: big compact blobs (bright stars) are down-weighted.

    Components smaller than ``area_cap`` (faint-trail fragments, faint stars)
    or thin ones keep weight 1; others get ``area_cap / area``, so a bright
    star cannot fake a line together with a few neighbours. "Thin" means
    ``area / max(bbox side)**2 < max_fill`` -- true for a trail, for crossing
    trails, for a trail touching a star or a bad column (where moment-based
    elongation fails), and false for round blobs (~0.79). Stars are not
    removed.
    """
    lab, n = label(hot, structure=np.ones((3, 3)))
    if n == 0:
        return np.zeros(0)
    yy, xx = np.nonzero(hot)
    lab_px = lab[yy, xx]
    area = np.bincount(lab_px, minlength=n + 1)[1:].astype(float)
    side = np.array([max(sl[0].stop - sl[0].start, sl[1].stop - sl[1].start)
                     for sl in find_objects(lab)], dtype=float)
    fill = area / side ** 2
    w = np.where((area <= area_cap) | (fill < max_fill), 1.0, area_cap / area)
    return w[lab_px - 1]


def _slab(p, u, lo, hi):
    """Parameter interval where p + t*u lies in [lo, hi] (vectorised)."""
    big = 1e12
    with np.errstate(divide="ignore", invalid="ignore"):
        t_a, t_b = (lo - p) / u, (hi - p) / u
    flat = np.abs(u) < 1e-12
    inside = (p >= lo) & (p <= hi)
    t_min = np.where(flat, np.where(inside, -big, big), np.minimum(t_a, t_b))
    t_max = np.where(flat, np.where(inside, big, -big), np.maximum(t_a, t_b))
    return t_min, t_max


def _chord(rho, theta, x0, x1, y0, y1):
    """(t_min, t_max) of lines inside the rectangle; t_max < t_min if missed."""
    c, s = np.cos(theta), np.sin(theta)
    tx_min, tx_max = _slab(rho * c, -s, x0, x1)
    ty_min, ty_max = _slab(rho * s, c, y0, y1)
    return np.maximum(tx_min, ty_min), np.minimum(tx_max, ty_max)


def _fold(rho, theta):
    if theta >= np.pi / 2:
        return -rho, theta - np.pi
    if theta < -np.pi / 2:
        return -rho, theta + np.pi
    return rho, theta


def _segment_signal(d, t, v, offset, hw, gap, flank, seg, t_lo, n_seg):
    """Per-segment signal of a line shifted by ``offset``.

    The band mean minus the flank median is taken on each side separately
    and the smaller one is returned: a trail is brighter than the background
    on *both* sides, a background step (e.g. an amplifier boundary) only on
    one.
    """
    dd = d - offset
    a = np.abs(dd)
    band = a <= hw
    side = (a > hw + gap) & (a <= hw + gap + flank)
    out = np.full(n_seg, np.nan)
    k_all = ((t - t_lo) // seg).astype(np.int64)
    ok = (k_all >= 0) & (k_all < n_seg)
    kb = k_all[band & ok]
    nb = np.bincount(kb, minlength=n_seg)
    sb = np.bincount(kb, weights=v[band & ok], minlength=n_seg)
    min_px = max(3, int(0.5 * seg * 2 * hw))
    min_side = max(3, min_px // 2)
    sides = []
    for sel in (side & ok & (dd < 0), side & ok & (dd > 0)):
        kf, vf = k_all[sel], v[sel]
        order = np.lexsort((vf, kf))
        kf, vf = kf[order], vf[order]
        bounds = np.searchsorted(kf, np.arange(n_seg + 1))
        med = np.full(n_seg, np.nan)
        for j in np.nonzero(np.diff(bounds) >= min_side)[0]:
            lo, hi = bounds[j], bounds[j + 1]
            m = (hi - lo) // 2
            med[j] = vf[lo + m] if (hi - lo) % 2 else 0.5 * (vf[lo + m - 1] + vf[lo + m])
        sides.append(med)
    good = (nb >= min_px) & np.isfinite(sides[0]) & np.isfinite(sides[1])
    mean_b = np.where(nb > 0, sb / np.maximum(nb, 1), np.nan)
    out[good] = np.minimum(mean_b[good] - sides[0][good], mean_b[good] - sides[1][good])
    return out


def _best_run(score):
    """(start, stop) of the maximum-sum contiguous run (Kadane)."""
    best, best_sum = None, 0.0
    cur, start = 0.0, 0
    for i, sc in enumerate(score):
        if cur <= 0:
            cur, start = sc, i
        else:
            cur += sc
        if cur > best_sum:
            best_sum, best = cur, (start, i + 1)
    return best


def _fwhm_profile(d, v, step=0.5, max_d=None):
    """FWHM and peak of a 1-D profile sampled at distances d (stacked pixels)."""
    if max_d is None:
        max_d = np.abs(d).max()
    edges = np.arange(-max_d, max_d + step, step)
    k = np.digitize(d, edges) - 1
    prof = np.array([np.median(v[k == i]) if np.any(k == i) else np.nan for i in range(len(edges) - 1)])
    centres = 0.5 * (edges[1:] + edges[:-1])
    ok = np.isfinite(prof)
    if ok.sum() < 5:
        return np.nan, np.nan
    prof, centres = prof[ok], centres[ok]
    i = int(np.argmax(prof))
    peak = prof[i]
    if peak <= 0:
        return np.nan, np.nan
    half = peak / 2
    left = i
    while left > 0 and prof[left] > half:
        left -= 1
    right = i
    while right < len(prof) - 1 and prof[right] > half:
        right += 1
    if prof[left] > half or prof[right] > half:
        return np.nan, peak

    def cross(a, b):
        y1, y2 = prof[a], prof[b]
        return centres[a] + (half - y1) / (y2 - y1) * (centres[b] - centres[a]) if y2 != y1 else centres[a]

    return cross(right - 1, right) - cross(left, left + 1), peak


# --- main entry point ------------------------------------------------------

def find_satellites(image, bin=None, k_input=2.0, z_min=8.0, nsigma=10.0,
                    min_fraction=0.3, min_length=100.0, width_guess=3.0,
                    bkg_box=64, edge=8, theta_step=0.25, max_candidates=300,
                    axis_tol=1.0, mask_nsigma=0.5, mask_margin=1.5, min_segments=10,
                    min_width=1.5, max_width=10.0, min_consistency=0.6, saturation=None,
                    pattern_min_lines=3, pattern_max_snr=2.0, min_peak_snr=0.5):
    """Find straight trails in ``image``.

    Parameters
    ----------
    image : 2-D array
    bin : int or None
        Binning factor for detection; ``None`` -> 1, 2 or 4 by frame size.
    k_input : float
        Threshold (sigma of the smoothed residual) for pixels entering Hough.
    z_min : float
        Minimal binomial significance of a Hough cell to become a candidate.
    nsigma : float
        Minimal significance of the median segment signal to accept a trail.
    min_fraction : float
        Minimal fraction of significant segments within the trail extent.
    min_consistency : float
        Minimal fraction of segments whose signal lies within
        max(3 sigma, 50%) of the median -- rejects lines crossing the halo
        of a bright star (a trail is roughly uniform along its length).
    min_length : float
        Minimal trail length in pixels.
    min_segments : int
        Minimal number of verification segments (~8 px each) in a trail;
        guards the median test against chains of a few stars.
    width_guess : float
        Expected trail FWHM in pixels; sets the verification band.
    bkg_box : int
        Background block size in pixels.
    edge : int
        Border (binned pixels) excluded from detection; trails reaching it
        are extended to the frame edge.
    theta_step : float
        Hough angle step in degrees (refined locally afterwards).
    max_candidates : int
        Maximal number of Hough candidates examined.
    axis_tol : float
        Lines within this many degrees of a detector axis are reported as
        "column"/"row" instead of "satellite".
    min_width, max_width : float
        Allowed FWHM (pixels) of a satellite trail. A trail is a smeared PSF:
        narrower features are detector structure or noise, much broader ones
        are low-level bands (flat or pattern structure). Columns/rows are
        exempt.
    saturation : float or None
        Saturation level (ADU). Saturated pixels (bleed trails, star cores)
        are excluded from voting and verification, so aligned bleeds of
        several stars are not taken for a trail.
    min_peak_snr : float
        Minimal trail peak in units of the pixel noise. Long, very faint
        linear structure (e.g. a flat-field pattern at 0.3-0.4 sigma) is
        statistically significant but is not a satellite; this is the
        sensitivity floor.
    pattern_min_lines, pattern_max_snr : int, float
        At least ``pattern_min_lines`` satellites parallel within 0.5 deg
        whose median ``peak_snr`` is below ``pattern_max_snr`` are a faint
        periodic pattern (flat/detector structure), reported as "pattern" and
        not masked. Bright parallel trails (e.g. a satellite train) stay
        satellites.
    mask_nsigma, mask_margin : float
        The mask covers the band where the trail profile exceeds
        ``mask_nsigma`` times the pixel noise (at least one FWHM from the
        centre), plus ``mask_margin`` pixels, also beyond the segment ends.

    Returns
    -------
    table : astropy.table.Table
        One row per detected linear feature, columns ``COLUMNS``.
    mask : 2-D bool array
        True on pixels covered by satellite trails (``kind == "satellite"``).
    """
    image = np.asarray(image)
    full_shape = image.shape
    b = _auto_bin(full_shape) if bin is None else int(bin)

    # --- 1. binned, background-flattened residual and the vote set ---------
    img = _bin(image, b)
    resid = img - _background(img, max(8, bkg_box // b))
    smooth = gaussian_filter(resid, 1.0)
    sig_s = mad_std(smooth)
    ny, nx = img.shape
    valid = np.zeros(img.shape, dtype=bool)
    valid[edge:ny - edge, edge:nx - edge] = True
    if saturation is not None:
        # saturated cores and bleed columns neither vote nor take part in
        # verification; a saturated trail is still seen through its wings
        sat_zone = binary_dilation(_bin(image >= saturation, b) > 0,
                                   iterations=max(1, int(round(4 / b))))
    else:
        sat_zone = np.zeros(img.shape, dtype=bool)
    hot = (smooth > k_input * sig_s) & valid & ~sat_zone
    ys, xs = np.nonzero(hot)
    ws = _component_weights(hot, area_cap=max(9.0, 60.0 / b ** 2))
    ys, xs = ys.astype(float), xs.astype(float)
    p = max(ws.sum() / max(valid.sum(), 1), 1e-6)

    # --- 2. Hough accumulator and its noise expectation --------------------
    thetas = np.deg2rad(np.arange(-90.0, 90.0, theta_step))
    cos_t, sin_t = np.cos(thetas), np.sin(thetas)
    rho_max = int(np.ceil(np.hypot(ny, nx)))
    rhos = np.arange(-rho_max, rho_max + 1, dtype=float)
    acc = _hough(ys, xs, ws, cos_t, sin_t, rho_max)
    box = (edge - 0.5, nx - edge - 0.5, edge - 0.5, ny - edge - 0.5)
    t_min, t_max = _chord(rhos[:, None], thetas[None, :], *box)
    chord = np.clip(t_max - t_min, 0, None)
    min_len_b = min_length / b
    usable = chord >= min_len_b
    expected = p * chord
    noise = np.sqrt(p * (1 - p) * chord + 1.0)

    # verification geometry (binned pixels)
    hw = max(1.5, width_guess / b)
    gap, flank = 2.0, 4.0
    seg = max(4.0, 8.0 / b)       # ~8 px (full resolution) per segment, at least 4 binned px
    step_off = hw + gap + flank + 3.0
    null_offsets = np.array([-1.0, 1.0, -2.0, 2.0]) * step_off
    reach = 2 * step_off + hw + gap + flank + 1.0

    # accumulator neighbourhood suppressed after a rejection: cells whose
    # lines share most pixels with the rejected one (a broad band spans many)
    dr_sup = int(np.ceil(hw + gap + 2))
    dt_sup = max(2, int(round(0.75 / theta_step)))

    yy, xx = np.indices(img.shape)
    xx, yy = xx.astype(float), yy.astype(float)
    open_px = valid.copy()        # pixels usable for verification
    open_px &= ~sat_zone

    results = []
    n_examined = 0
    log = []
    z = np.where(usable, (acc - expected) / noise, -np.inf)
    for _ in range(max_candidates):
        n_examined += 1
        ri, ti = np.unravel_index(int(np.argmax(z)), z.shape)
        if z[ri, ti] < z_min:
            break
        rho, theta = rhos[ri], thetas[ti]
        cand = dict(z=float(z[ri, ti]), rho_b=float(rho), theta_deg=float(np.rad2deg(theta)),
                    outcome="rejected")
        log.append(cand)

        # --- 3a. local refinement on the vote set ---------------------------
        d_v = xs * np.cos(theta) + ys * np.sin(theta) - rho
        near = np.abs(d_v) <= 4.0
        if near.sum() >= 3:
            xn, yn = xs[near], ys[near]
            best = (-1, rho, theta)
            for dth in np.deg2rad(np.linspace(-theta_step, theta_step, 13)):
                th = theta + dth
                dd = xn * np.cos(th) + yn * np.sin(th) - rho
                cnt = np.array([np.sum(np.abs(dd - dr) <= 1.0) for dr in np.arange(-2.0, 2.01, 0.25)])
                j = int(np.argmax(cnt))
                if cnt[j] > best[0]:
                    best = (cnt[j], rho - 2.0 + 0.25 * j, th)
            _, rho, theta = best

        # --- 3b. verification on the residual image -------------------------
        c, s = np.cos(theta), np.sin(theta)
        d_full = xx * c + yy * s - rho
        sub = open_px & (np.abs(d_full) <= reach)
        accepted = False
        if sub.sum() > 0:
            d = d_full[sub]
            t = -xx[sub] * s + yy[sub] * c
            v = resid[sub]
            t_lo, t_hi = _chord(rho, theta, *box)
            n_seg = int(np.ceil((t_hi - t_lo) / seg)) if t_hi > t_lo else 0
            if n_seg >= 3:
                sig_v = _segment_signal(d, t, v, 0.0, hw, gap, flank, seg, t_lo, n_seg)
                nulls = np.concatenate([_segment_signal(d, t, v, off, hw, gap, flank, seg, t_lo, n_seg)
                                        for off in null_offsets])
                nulls = nulls[np.isfinite(nulls)]
                sig_seg = mad_std(nulls) if nulls.size >= 10 else np.nan
                if np.isfinite(sig_seg) and sig_seg > 0:
                    s_k = sig_v / sig_seg
                    score = np.where(np.isfinite(s_k), np.clip(s_k, -4.0, 4.0) - 1.0, 0.0)
                    run = _best_run(score)
                    if run is not None:
                        i0, i1 = run
                        sv = sig_v[i0:i1]
                        sv = sv[np.isfinite(sv)]
                        if sv.size >= min_segments and (i1 - i0) * seg >= min_len_b:
                            med = np.median(sv)
                            significance = med / (1.2533 * sig_seg / np.sqrt(sv.size))
                            fraction = np.mean(sv > 2.0 * sig_seg)
                            # a trail is roughly uniform along its length; a
                            # line through a bright star's halo peaks mid-way
                            consistency = np.mean(np.abs(sv - med) <= max(3.0 * sig_seg, 0.5 * abs(med)))
                            cand.update(length_b=float((i1 - i0) * seg), significance=float(significance),
                                        fraction=float(fraction), consistency=float(consistency))
                            if (significance >= nsigma and fraction >= min_fraction
                                    and consistency >= min_consistency):
                                accepted = True
                                finite = np.nonzero(np.isfinite(sig_v))[0]
                                ext = (t_lo + i0 * seg, t_lo + i1 * seg)
                                to_edge = (i0 <= finite[0] + 1, i1 >= finite[-1])

        if not accepted:
            # suppress this neighbourhood of the accumulator and move on
            z[max(0, ri - dr_sup):ri + dr_sup + 1, max(0, ti - dt_sup):ti + dt_sup + 1] = -np.inf
            continue

        # --- 4. full-resolution refinement ---------------------------------
        row = _refine_full(image, b, rho, theta, ext, to_edge, hw, seg, mask_nsigma, mask_margin)
        row["significance"] = float(significance)
        row["fraction"] = float(fraction)
        deg = np.rad2deg(row["theta"])
        if abs(deg) < axis_tol:
            row["kind"] = "column"
        elif abs(abs(deg) - 90.0) < axis_tol:
            row["kind"] = "row"
        elif (np.isfinite(row["peak"]) and min_width <= row["width"] <= max_width
              and row["peak_snr"] >= min_peak_snr):
            row["kind"] = "satellite"
        else:
            # narrower than any PSF-smeared trail, a broad band, or too faint
            # to tell from low-level flat/pattern structure: not a satellite
            cand["outcome"] = "shape"
            z[max(0, ri - dr_sup):ri + dr_sup + 1, max(0, ti - dt_sup):ti + dt_sup + 1] = -np.inf
            continue
        cand["outcome"] = row["kind"]
        results.append(row)

        # remove this feature's votes and pixels, and continue
        t_px = -xx * s + yy * c
        open_px &= ~((np.abs(d_full) <= hw + gap) & (t_px >= ext[0] - seg) & (t_px <= ext[1] + seg))
        d_v = xs * c + ys * s - rho
        t_v = -xs * s + ys * c
        gone = (np.abs(d_v) <= hw + 2.0) & (t_v >= ext[0] - seg) & (t_v <= ext[1] + seg)
        if gone.any():
            acc -= _hough(ys[gone], xs[gone], ws[gone], cos_t, sin_t, rho_max)
            ys, xs, ws = ys[~gone], xs[~gone], ws[~gone]
            suppressed = ~np.isfinite(z)
            z = np.where(usable & ~suppressed, (acc - expected) / noise, -np.inf)

    _mark_patterns(results, pattern_min_lines, pattern_max_snr)
    table = Table(rows=[[r[c] for c in COLUMNS] for r in results] if results else None,
                  names=COLUMNS, dtype=[float] * (len(COLUMNS) - 1) + [str])
    # bookkeeping for tuning: did the candidate budget run out?
    table.meta.update(bin=b, candidates_examined=n_examined,
                      budget_exhausted=bool(n_examined >= max_candidates), candidates=log)
    mask = np.zeros(full_shape, dtype=bool)
    for r in results:
        if r["kind"] == "satellite":
            _paint_segment(mask, r, mask_margin)
    return table, mask


def _group_medians(k, v, n):
    """Median of v for each group index k in [0, n); NaN where empty."""
    out = np.full(n, np.nan)
    if k.size == 0:
        return out
    order = np.lexsort((v, k))
    k, v = k[order], v[order]
    bounds = np.searchsorted(k, np.arange(n + 1))
    for j in np.nonzero(np.diff(bounds) > 0)[0]:
        lo, hi = bounds[j], bounds[j + 1]
        m = (hi - lo) // 2
        out[j] = v[lo + m] if (hi - lo) % 2 else 0.5 * (v[lo + m - 1] + v[lo + m])
    return out


def _refine_full(image, b, rho_b, theta, ext_b, to_edge, hw_b, seg_b, mask_nsigma, mask_margin):
    """Refine a binned detection on the full-resolution image."""
    ny, nx = image.shape
    c, s = np.cos(theta), np.sin(theta)
    off = (b - 1) / 2.0                       # binned pixel centre in full-res px
    rho = b * rho_b + off * (c + s)
    t0 = b * ext_b[0] + off * (c - s)
    t1 = b * ext_b[1] + off * (c - s)
    seg = max(8.0, seg_b * b)
    search = max(6.0, (hw_b + 2) * b)
    outer = search + 8.0
    frame = (-0.5, nx - 0.5, -0.5, ny - 0.5)

    def ends(rho, theta, t0, t1):
        c, s = np.cos(theta), np.sin(theta)
        return np.array([[rho * c - t * s, rho * s + t * c] for t in (t0, t1)])

    if to_edge[0] or to_edge[1]:
        f0, f1 = _chord(rho, theta, *frame)
        t0 = f0 if to_edge[0] else t0
        t1 = f1 if to_edge[1] else t1

    # pixels near the segment (bounding box first, then a distance cut)
    pts = ends(rho, theta, t0, t1)
    pad = outer + 4
    x_lo, x_hi = int(max(0, pts[:, 0].min() - pad)), int(min(nx, pts[:, 0].max() + pad + 1))
    y_lo, y_hi = int(max(0, pts[:, 1].min() - pad)), int(min(ny, pts[:, 1].max() + pad + 1))
    yy, xx = np.mgrid[y_lo:y_hi, x_lo:x_hi]
    d0 = xx * c + yy * s - rho
    near = np.abs(d0) <= outer + 4
    X = xx[near].astype(float)
    Y = yy[near].astype(float)
    V = image[y_lo:y_hi, x_lo:x_hi][near].astype(float)

    def geometry(rho, theta, t0, t1, window=search):
        c, s = np.cos(theta), np.sin(theta)
        d = X * c + Y * s - rho
        t = -X * s + Y * c
        inb = (t >= t0) & (t <= t1)
        core = inb & (np.abs(d) <= window)
        bkgz = inb & (np.abs(d) > search + 2) & (np.abs(d) <= outer)
        return d, t, core, bkgz

    # Matched filter: maximise the mean (clipped) signal in a narrow band
    # over a small (dtheta, drho) grid -- uses the whole length, so it is
    # accurate also for faint trails; clipping keeps stars from deciding.
    d, t, core, bkgz = geometry(rho, theta, t0, t1)
    if bkgz.sum() >= 10:
        bkg, noise = np.median(V[bkgz]), mad_std(V[bkgz])
        inb = (t >= t0) & (t <= t1) & (np.abs(d) <= search)
        Vc = np.clip(V[inb] - bkg, -5 * noise, 5 * noise)
        Xi, Yi = X[inb], Y[inb]
        length = max(t1 - t0, 1.0)
        dth_max = min(np.deg2rad(0.5), 2.0 * b / length)
        step = 0.25
        edges = np.arange(-search, search + step, step)
        nb = max(1, int(round(1.5 / step)))          # +-1.5 px band = 12 bins
        best = (-np.inf, 0.0, 0.0)
        for dth in np.linspace(-dth_max, dth_max, 21):
            th = theta + dth
            dd = Xi * np.cos(th) + Yi * np.sin(th) - rho
            sums, _ = np.histogram(dd, edges, weights=Vc)
            cnts, _ = np.histogram(dd, edges)
            cs, cc = np.r_[0, np.cumsum(sums)], np.r_[0, np.cumsum(cnts)]
            band_s = cs[2 * nb:] - cs[:-2 * nb]
            band_c = cc[2 * nb:] - cc[:-2 * nb]
            score = np.where(band_c > 0, band_s / np.maximum(band_c, 1), -np.inf)
            j = int(np.argmax(score))
            if score[j] > best[0]:
                best = (score[j], dth, edges[j + nb])
        _, dth, dr = best
        pts = ends(rho, theta, t0, t1)
        theta, rho = theta + dth, rho + dr
        c, s = np.cos(theta), np.sin(theta)
        t_ends = pts[:, 0] * -s + pts[:, 1] * c
        t0, t1 = t_ends.min(), t_ends.max()

    window = max(2.0, 0.75 * hw_b * b)
    for _ in range(2):                         # centroid refinement, narrow window
        d, t, core, bkgz = geometry(rho, theta, t0, t1, window)
        n_seg = max(1, int((t1 - t0) // seg))
        k = np.minimum(((t - t0) // seg).astype(np.int64), n_seg - 1)
        bkg = _group_medians(k[bkgz], V[bkgz], n_seg)
        kc = k[core]
        w = np.clip(V[core] - bkg[kc], 0, None)
        w = np.where(np.isfinite(w), w, 0.0)
        sw = np.bincount(kc, weights=w, minlength=n_seg)
        swd = np.bincount(kc, weights=w * d[core], minlength=n_seg)
        swt = np.bincount(kc, weights=w * t[core], minlength=n_seg)
        good = sw > 0
        if good.sum() < 3:
            break
        tc, dc = swt[good] / sw[good], swd[good] / sw[good]
        beta, a = np.polyfit(tc, dc, 1, w=np.sqrt(sw[good]))
        # the new line is d - (a + beta t) = 0; keep the same physical ends
        pts = ends(rho, theta, t0, t1)
        theta = theta - np.arctan(beta)
        rho = (rho + a) / np.hypot(1.0, beta)
        c, s = np.cos(theta), np.sin(theta)
        t_ends = pts[:, 0] * -s + pts[:, 1] * c
        t0, t1 = t_ends.min(), t_ends.max()
        if to_edge[0] or to_edge[1]:
            f0, f1 = _chord(rho, theta, *frame)
            t0 = f0 if to_edge[0] else t0
            t1 = f1 if to_edge[1] else t1

    # fold theta into [-pi/2, pi/2); t flips sign with it, so carry the
    # physical end points across
    pts = ends(rho, theta, t0, t1)
    rho, theta = _fold(rho, theta)
    c, s = np.cos(theta), np.sin(theta)
    t_ends = pts[:, 0] * -s + pts[:, 1] * c
    t0, t1 = t_ends.min(), t_ends.max()

    d, t, core, bkgz = geometry(rho, theta, t0, t1)
    if bkgz.sum() >= 10:
        bkg, noise = np.median(V[bkgz]), mad_std(V[bkgz])
    else:
        bkg, noise = np.median(V), mad_std(V)
    # width/peak stay NaN when the profile cannot be measured; such a
    # feature cannot be confirmed as a trail (columns/rows are exempt)
    width, peak = _fwhm_profile(d[core], V[core] - bkg, step=0.5, max_d=search)
    w_eff = width if np.isfinite(width) else float(hw_b * b)
    half = w_eff
    if np.isfinite(peak) and noise > 0 and peak > mask_nsigma * noise:
        sigma_w = w_eff / 2.3548
        half = max(w_eff, sigma_w * np.sqrt(2.0 * np.log(peak / (mask_nsigma * noise))))
    peak_snr = peak / noise if (np.isfinite(peak) and noise > 0) else np.nan
    (x0, y0), (x1, y1) = ends(rho, theta, t0, t1)
    return dict(x0=float(x0), y0=float(y0), x1=float(x1), y1=float(y1),
                rho=float(rho), theta=float(theta), length=float(t1 - t0),
                width=float(width), peak=float(peak) if np.isfinite(peak) else np.nan,
                peak_snr=float(peak_snr), mask_halfwidth=float(half + mask_margin))


def _mark_patterns(results, min_lines, max_snr, tol_deg=0.5):
    """Relabel groups of faint parallel 'satellites' as 'pattern'."""
    sats = [r for r in results if r["kind"] == "satellite"]
    used = set()
    for i, r in enumerate(sats):
        if i in used:
            continue
        group = [j for j, q in enumerate(sats)
                 if np.rad2deg(abs((q["theta"] - r["theta"] + np.pi / 2) % np.pi - np.pi / 2)) < tol_deg]
        if len(group) >= min_lines and np.nanmedian([sats[j]["peak_snr"] for j in group]) < max_snr:
            for j in group:
                sats[j]["kind"] = "pattern"
            used.update(group)


def _paint_segment(mask, row, mask_margin):
    ny, nx = mask.shape
    c, s = np.cos(row["theta"]), np.sin(row["theta"])
    half = row["mask_halfwidth"]
    t0 = -row["x0"] * s + row["y0"] * c
    t1 = -row["x1"] * s + row["y1"] * c
    t0, t1 = min(t0, t1) - mask_margin, max(t0, t1) + mask_margin
    pad = half + mask_margin + 1
    x_lo = int(max(0, min(row["x0"], row["x1"]) - pad))
    x_hi = int(min(nx, max(row["x0"], row["x1"]) + pad + 1))
    y_lo = int(max(0, min(row["y0"], row["y1"]) - pad))
    y_hi = int(min(ny, max(row["y0"], row["y1"]) + pad + 1))
    yy, xx = np.mgrid[y_lo:y_hi, x_lo:x_hi]
    d = xx * c + yy * s - row["rho"]
    t = -xx * s + yy * c
    mask[y_lo:y_hi, x_lo:x_hi] |= (np.abs(d) <= half) & (t >= t0) & (t <= t1)
