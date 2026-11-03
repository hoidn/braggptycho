#!/usr/bin/env python3
"""First-look analysis for 26-ID nanoprobe Bragg scan files (NbN thin films).

Input:  scan h5 per David's email: /entry/data/data = stack of 2D detector
        frames (nth frame <-> nth ZPX value), /entry/scan/Positioner/{ZPX,ZPY}.
Output: <prefix>_maps.npz + <prefix>_quicklook.png with
        - mean diffraction pattern (log) and an example frame
        - virtual-detector maps: total intensity, max pixel
        - Bragg-peak COM maps (detector x/y, per-frame, ROI-based)
        - adjacent-pattern Pearson correlation maps (fast/slow axis)
        - correlation decay vs step offset along both axes ("4 numbers" idea)

Usage:  python quicklook.py scan_308_100308.h5 [--roi-half 64] [--out prefix]
        python quicklook.py --selftest
"""
import argparse
import sys

import numpy as np

try:
    import hdf5plugin  # noqa: F401  (registers LZ4/bslz4 codecs used by detector files)
except ImportError:
    pass
import h5py
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle

FLAG = 4294967295  # 2**32-1: detector bad-pixel / module-gap sentinel in 26-ID uint32 frames


def find_frames(f):
    """Return the frame stack dataset (lazy h5py dataset, shape (N, ny, nx))."""
    if "entry/data/data" in f:
        return f["entry/data/data"]
    hits = []
    f.visititems(lambda n, o: hits.append(o) if isinstance(o, h5py.Dataset) and o.ndim == 3 else None)
    if not hits:
        raise KeyError("no 3D frame stack found in file")
    return max(hits, key=lambda d: d.size)


def find_positions(f):
    """Return (zpx, zpy) flattened to per-frame arrays if present, else (None, None)."""
    zp = {}
    def visit(name, obj):
        if isinstance(obj, h5py.Dataset) and name.rsplit("/", 1)[-1].upper() in ("ZPX", "ZPY"):
            zp[name.rsplit("/", 1)[-1].upper()] = obj[...]
    f.visititems(visit)
    return zp.get("ZPX"), zp.get("ZPY")


def _factor_pairs(n):
    """(rows, cols) factor pairs, most-square first."""
    return [(k, n // k) for k in range(int(np.sqrt(n)), 0, -1) if n % k == 0]


def infer_grid(n_frames, zpx, zpy, grid=None):
    """Infer (n_rows, n_cols); ZPX is the fast (inner) axis.

    Priority: explicit grid > per-frame positioner arrays > factorization
    (orientation resolved later by map coherence in run()).
    """
    if grid is not None:
        assert grid[0] * grid[1] == n_frames, f"grid {grid} != {n_frames} frames"
        return tuple(grid)
    if zpy is not None and zpy.size == n_frames:
        n_rows = len(np.unique(np.round(zpy.ravel(), 6)))
        if n_rows > 1 and n_frames % n_rows == 0:
            return n_rows, n_frames // n_rows
    pairs = _factor_pairs(n_frames)
    if len(pairs) == 1 and pairs[0][0] == 1:
        raise ValueError(f"{n_frames} is prime; pass --grid R,C")
    return pairs[0]  # provisional; run() may transpose by coherence


def process(frames, n_rows, n_cols, roi_half):
    """Single pass over frames -> stats + ROI crops (float32, in memory).

    Flagged pixels (== FLAG) are masked to 0 so they don't corrupt the peak
    search, per-frame sums, or COM.
    """
    n, fy, fx = frames.shape
    step = max(1, n // 80)
    acc = np.zeros((fy, fx), dtype=np.float64)
    cnt = np.zeros((fy, fx), dtype=np.float64)
    for i in range(0, n, step):
        fr = frames[i]
        good = fr < FLAG
        acc[good] += fr[good]
        cnt[good] += 1
    mean_pat = np.where(cnt > 0, acc / np.maximum(cnt, 1), 0.0)
    cy, cx = np.unravel_index(np.argmax(mean_pat), mean_pat.shape)
    y0, y1 = max(0, cy - roi_half), min(fy, cy + roi_half)
    x0, x1 = max(0, cx - roi_half), min(fx, cx + roi_half)
    ry, rx = y1 - y0, x1 - x0
    est_gb = n * ry * rx * 4 / 1e9
    if est_gb > 8:
        raise MemoryError(f"ROI stack would need {est_gb:.1f} GB; shrink --roi-half")

    crops = np.empty((n, ry, rx), dtype=np.float32)
    total = np.empty(n)
    peak = np.empty(n)
    comy = np.empty(n)
    comx = np.empty(n)
    iy, ix = np.mgrid[0:ry, 0:rx]
    chunk = max(1, min(n, int(64e6 // (fy * fx * 4))))
    roi_good = None
    for s in range(0, n, chunk):
        block = frames[s:s + chunk]
        c = np.asarray(block[:, y0:y1, x0:x1], dtype=np.float32)
        if roi_good is None:
            roi_good = c[0] < FLAG  # bad-pixel/module-gap layout is static
        c[c >= FLAG] = 0.0  # mask bad pixels inside the ROI
        crops[s:s + chunk] = c
        t = c.sum(axis=(1, 2))
        total[s:s + chunk] = t
        peak[s:s + chunk] = c.max(axis=(1, 2))
        safe = np.where(t > 0, t, 1.0)
        comy[s:s + chunk] = (c * iy).sum(axis=(1, 2)) / safe - (ry - 1) / 2
        comx[s:s + chunk] = (c * ix).sum(axis=(1, 2)) / safe - (rx - 1) / 2

    shape = (n_rows, n_cols)
    return dict(
        mean_pattern=mean_pat, roi=(y0, y1, x0, x1), crops=crops, roi_good=roi_good,
        total=total.reshape(shape), peak=peak.reshape(shape),
        comy=comy.reshape(shape), comx=comx.reshape(shape),
    )


def pearson_maps(crops, n_rows, n_cols, max_offset=8, good=None):
    """Adjacent-pattern Pearson correlation along fast (x) and slow (y) axes.

    Returns a dict of maps (NaN on last col/row) and decay curves (offset
    1..max_offset), in raw and Poisson-noise-corrected form. For photon-count
    frames shot noise adds variance equal to the per-frame mean, so
    r_corrected = cov(a,b)/sqrt((var a - mean a)(var b - mean b)) estimates the
    correlation of the underlying signals (unbounded within noise; clamp on
    display, not here). `good` restricts to unflagged detector pixels.
    """
    n = crops.shape[0]
    flat = crops.reshape(n, -1)
    if good is not None:
        flat = flat[:, np.asarray(good).ravel()]
    g = flat.astype(np.float64)
    mean = g.mean(axis=1)
    var = g.var(axis=1)
    sigvar = np.maximum(var - mean, 1e-9)  # Poisson: noise variance == mean counts
    g -= mean[:, None]
    idx = np.arange(n).reshape(n_rows, n_cols)

    def pair(i, j):
        cov = (g[i] * g[j]).mean(axis=-1)
        return cov / np.sqrt(var[i] * var[j]), cov / np.sqrt(sigvar[i] * sigvar[j])

    out = {k: np.full((n_rows, n_cols), np.nan) for k in ("corr_x", "corr_y", "corrc_x", "corrc_y")}
    out["corr_x"][:, :-1], out["corrc_x"][:, :-1] = (
        m.reshape(n_rows, n_cols - 1) for m in pair(idx[:, :-1].ravel(), idx[:, 1:].ravel()))
    out["corr_y"][:-1, :], out["corrc_y"][:-1, :] = (
        m.reshape(n_rows - 1, n_cols) for m in pair(idx[:-1].ravel(), idx[1:].ravel()))
    for axis, tag in ((1, "x"), (0, "y")):
        raw, cor = [], []
        for k in range(1, max_offset + 1):
            i = idx.take(range(idx.shape[axis] - k), axis=axis).ravel()
            j = idx.take(range(k, idx.shape[axis]), axis=axis).ravel()
            r, c = pair(i, j)
            raw.append(np.median(r))
            cor.append(np.median(c))
        out["decay_" + tag] = np.array(raw)
        out["decayc_" + tag] = np.array(cor)
    return out


def _roughness(a):
    """Mean abs adjacent difference of a z-scored map; lower = more coherent."""
    a = (a - a.mean()) / (a.std() + 1e-9)
    return float(np.mean(np.abs(np.diff(a, axis=1))) + np.mean(np.abs(np.diff(a, axis=0))))


def run(path, roi_half=64, out=None, max_offset=8, grid=None):
    out = out or path.rsplit(".", 1)[0]
    with h5py.File(path, "r") as f:
        frames = find_frames(f)
        zpx, zpy = find_positions(f)
        n_rows, n_cols = infer_grid(frames.shape[0], zpx, zpy, grid)
        print(f"{path}: {frames.shape[0]} frames {frames.shape[1]}x{frames.shape[2]}, "
              f"grid {n_rows}x{n_cols} (provisional)")
        r = process(frames, n_rows, n_cols, roi_half)
    # resolve square-ambiguous orientation by which reshape gives a smoother intensity map
    if grid is None and n_rows != n_cols:
        alt = r["total"].reshape(n_cols, n_rows)
        if _roughness(alt) < _roughness(r["total"]):
            n_rows, n_cols = n_cols, n_rows
            for k in ("total", "peak", "comy", "comx"):
                r[k] = r[k].reshape(n_rows, n_cols)
            print(f"  -> transposed to grid {n_rows}x{n_cols} (more coherent)")
    pm = pearson_maps(r["crops"], n_rows, n_cols, max_offset, good=r["roi_good"])

    np.savez_compressed(
        out + "_maps.npz", mean_pattern=r["mean_pattern"], roi=np.array(r["roi"]),
        total=r["total"], peak=r["peak"], comx=r["comx"], comy=r["comy"],
        zpx=np.asarray(zpx) if zpx is not None else np.array([]),
        zpy=np.asarray(zpy) if zpy is not None else np.array([]),
        **pm,
    )

    fig, ax = plt.subplots(2, 4, figsize=(17, 8))
    y0, y1, x0, x1 = r["roi"]
    ax[0, 0].imshow(np.log1p(r["mean_pattern"]), cmap="magma")
    ax[0, 0].add_patch(Rectangle((x0, y0), x1 - x0, y1 - y0, ec="cyan", fc="none"))
    ax[0, 0].set_title("mean pattern (log), ROI")
    mid = r["crops"][len(r["crops"]) // 2]
    ax[0, 1].imshow(np.log1p(mid), cmap="magma")
    ax[0, 1].set_title("example frame (log, ROI)")
    for a, key, cmap in ((ax[0, 2], "total", "viridis"), (ax[0, 3], "peak", "viridis"),
                         (ax[1, 0], "comx", "RdBu_r"), (ax[1, 1], "comy", "RdBu_r")):
        im = a.imshow(r[key], cmap=cmap)
        a.set_title(key + " map")
        fig.colorbar(im, ax=a, shrink=0.8)
    for a, m, name in ((ax[1, 2], pm["corr_y"], "corr slow(y) raw"),
                       (ax[1, 3], pm["corrc_y"], "corr slow(y) noise-corrected")):
        im = a.imshow(m, cmap="inferno", vmin=0.5, vmax=1)
        a.set_title(name)
        fig.colorbar(im, ax=a, shrink=0.8)
    fig.suptitle(f"{path}  |  grid {r['total'].shape[0]}x{r['total'].shape[1]}  |  "
                 f"raw decay_x={np.round(pm['decay_x'][:4], 3)}  decay_y={np.round(pm['decay_y'][:4], 3)}  |  "
                 f"corrected x0={pm['decayc_x'][0]:.3f} y0={pm['decayc_y'][0]:.3f}")
    fig.tight_layout()
    fig.savefig(out + "_quicklook.png", dpi=130)
    print("wrote", out + "_maps.npz", "and", out + "_quicklook.png")
    return r, pm


def selftest(path="synthetic_scan.h5"):
    """Generate a synthetic twin-domain Bragg scan in David's h5 layout and verify."""
    rng = np.random.default_rng(0)
    npx, pix_nm, probe_fwhm_nm, step_nm, grid = 512, 10.0, 150.0, 30.0, 41
    # twin domains ~300 nm (>> illumination) so a scan position sits mostly in one variant
    blobs = np.fft.ifft2(np.fft.fft2(rng.standard_normal((npx, npx))) *
                         np.exp(-np.hypot(*np.mgrid[-npx//2:npx//2, -npx//2:npx//2]) ** 2
                                / (2 * 30.0 ** 2))).real
    domain = blobs > 0
    amp = np.where(domain, 1.0, 0.35)  # (202) lights variant A; B is a weak amplitude floor
    # fine-scale disorder (~30 nm correlation length) -> genuine speckle that decorrelates as the probe steps
    fine = np.fft.ifft2(np.fft.fft2(rng.standard_normal((npx, npx))) *
                        np.exp(-np.hypot(*np.mgrid[-npx//2:npx//2, -npx//2:npx//2]) ** 2
                               / (2 * 3.0 ** 2))).real
    obj = amp * np.exp(1j * 1.6 * fine / fine.std())

    fov = 96
    sig = probe_fwhm_nm / pix_nm / 2.355
    py, px = np.mgrid[-fov//2:fov//2, -fov//2:fov//2]
    probe = np.exp(-(py ** 2 + px ** 2) / (2 * sig ** 2))
    # domain-correlated bulk lattice tilt -> rigid shift of the Bragg speckle (mosaic = COM shift)
    tiltA, tiltB = 0.10, -0.10  # rad/pixel of local phase slope
    lx = px - px.mean()

    step_px = step_nm / pix_nm
    centers = npx // 2 + (np.arange(grid) - grid // 2) * step_px
    frames, zpx_list, zpy_list = [], [], []
    for r_i, cy in enumerate(centers):
        for c_i, cx in enumerate(centers):
            iy, ix = int(round(cy)), int(round(cx))
            patch = obj[iy - fov//2:iy + fov//2, ix - fov//2:ix + fov//2]
            slope = tiltA if domain[iy, ix] else tiltB
            wave = probe * patch * np.exp(1j * slope * lx)
            far = np.abs(np.fft.fftshift(np.fft.fft2(wave))) ** 2
            far *= 3e4 / far.sum()  # high flux: this is a code self-test, not the photon-starved regime
            frames.append(rng.poisson(far).astype(np.uint16))
            zpx_list.append((c_i - grid // 2) * step_nm / 1000)
            zpy_list.append((r_i - grid // 2) * step_nm / 1000)

    with h5py.File(path, "w") as f:
        f.create_dataset("entry/data/data", data=np.array(frames), chunks=(1, fov, fov),
                         compression="gzip", compression_opts=1)
        f.create_dataset("entry/scan/Positioner/ZPX", data=np.array(zpx_list))
        f.create_dataset("entry/scan/Positioner/ZPY", data=np.array(zpy_list))

    r, pm = run(path, roi_half=48)
    decay_x, decay_y = pm["decay_x"], pm["decay_y"]
    assert r["total"].shape == (grid, grid)
    assert np.isfinite(decay_x).all() and np.isfinite(decay_y).all()
    # aggregate decay (near pairs more correlated than far); element-wise monotonicity is not
    # guaranteed because the per-domain tilt flip makes boundary-straddling offset-1 pairs dip
    assert decay_x[:2].mean() > decay_x[-2:].mean(), f"corr should decay: {decay_x.round(3)}"
    # shot noise attenuates raw correlation; the Poisson correction must recover some of it
    assert pm["decayc_x"][0] > decay_x[0], f"corrected {pm['decayc_x'][0]:.3f} !> raw {decay_x[0]:.3f}"
    dom_map = domain[np.ix_(np.round(centers).astype(int), np.round(centers).astype(int))]
    com_split = abs(np.nanmean(r["comx"][dom_map]) - np.nanmean(r["comx"][~dom_map]))
    assert com_split > 1.0, f"twin-domain COM contrast not detected ({com_split:.2f} px)"
    print(f"SELFTEST PASS: grid ok, decay {decay_x[:2].mean():.3f}->{decay_x[-2:].mean():.3f}, "
          f"corrected corr {pm['decayc_x'][0]:.3f} > raw {decay_x[0]:.3f}, "
          f"domain COM split {com_split:.1f} px")


if __name__ == "__main__":
    p = argparse.ArgumentParser()
    p.add_argument("h5file", nargs="?", help="scan .h5 file")
    p.add_argument("--roi-half", type=int, default=64)
    p.add_argument("--out", default=None)
    p.add_argument("--grid", default=None, help="rows,cols (e.g. 41,71) to skip inference")
    p.add_argument("--selftest", action="store_true")
    a = p.parse_args()
    grid = tuple(int(x) for x in a.grid.split(",")) if a.grid else None
    if a.selftest:
        selftest()
    elif a.h5file:
        run(a.h5file, a.roi_half, a.out, grid=grid)
    else:
        p.print_help()
        sys.exit(1)
