#!/usr/bin/env python3
"""The 2D cross-correlation of neighbouring speckle IMAGES, as a function of lag.

Not the Bragg peak's rigid displacement -- the correlation surface of the
speckle itself:

    C_k(dq) = < r_i(q) r_j(q + dq) >     for scan neighbours j = i + k steps

with r = frame minus the flux-matched scan-mean pattern, so the smooth peak
envelope is removed and only the interference structure enters. k = 0 is the
autocorrelation (speckle size and shape); k >= 1 are the cross-correlations
(how that structure survives, deforms and decays as the probe steps).

This is the full 2D object the meeting's "correlate the two images next to each
other" asks for. Reading it:
  - PEAK HEIGHT vs k  -> decorrelation, i.e. how fast new material enters the
    coherent volume (the length scales quantified in speckle_geometry.py);
  - WIDTH/SHAPE       -> the coherently illuminated extent, and whether it
    changes with k (it should not, if the probe is stable);
  - CENTROID          -> any systematic translation of the speckle;
  - the zero-lag SPIKE is white shot noise and appears ONLY at k = 0 (noise is
    independent between exposures), which is why the k >= 1 surfaces are clean
    and give a noise-free view of the speckle shape.

Usage: python xcorr_2d.py [--selftest]
"""
import argparse

import numpy as np
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt

from quicklook import find_frames, process
from speckle_geometry import SCANS, extent_nm

try:
    import hdf5plugin  # noqa: F401  (registers bslz4)
except ImportError:
    pass
import h5py

N_ROWS, N_COLS = 71, 41
STEP_NM = 20.0
HALF = 14                 # lag window shown, px


def residual_ffts(crops, good):
    """Padded rfft2 of each envelope-subtracted frame and its corrected signal power.

    S_i = sum(r_i^2) - sum(counts_i) over good pixels: the shot-noise contribution to
    the residual power is the total counts (Poisson variance == mean), so removing it
    leaves the true speckle power. Normalising the cross-correlation by sqrt(S_i S_j)
    turns the surface peak into an honest correlation COEFFICIENT -- 1.0 at zero lag
    for a frame with itself, and directly comparable between the two reflections
    despite their ~17x difference in raw counts.
    """
    mean_pat = crops.mean(axis=0) * good
    denom = mean_pat.sum()
    ny, nx = mean_pat.shape
    out = np.empty((len(crops), 2 * ny, nx + 1), dtype=np.complex64)
    S = np.empty(len(crops))
    for i, f in enumerate(crops):
        c = f * good
        r = c - (c.sum() / denom) * mean_pat
        out[i] = np.fft.rfft2(r, s=(2 * ny, 2 * nx)).astype(np.complex64)
        S[i] = max((r ** 2).sum() - c.sum(), 1e-9)
    return out, (ny, nx), S


def mean_xcorr(F, shape, pairs, S=None, half=HALF, good=None):
    """Mean linear 2D cross-correlation surface over (i, j) frame pairs.

    With S supplied each pair is normalised by sqrt(S_i S_j) before averaging, so the
    surface is a correlation coefficient rather than an arbitrary power. The padded FFT
    prevents edge wraparound; mask-overlap normalisation removes the finite-window taper.
    """
    ny, nx = shape
    fft_shape = (2 * ny, 2 * nx)
    cross_power = np.zeros(F.shape[1:], dtype=np.complex128)
    for i, j in pairs:
        c = np.conj(F[i]) * F[j]
        cross_power += c / np.sqrt(S[i] * S[j]) if S is not None else c
    acc = np.fft.fftshift(np.fft.irfft2(cross_power / len(pairs), s=fft_shape))
    window = np.ones(shape) if good is None else np.asarray(good, dtype=float)
    W = np.fft.rfft2(window, s=fft_shape)
    overlap = np.fft.fftshift(np.fft.irfft2(np.abs(W) ** 2, s=fft_shape))
    acc *= window.sum() / np.maximum(overlap, 1)
    cy, cx = ny, nx
    return acc[cy - half:cy + half + 1, cx - half:cx + half + 1]


def neighbour_pairs(k, axis, n_rows=N_ROWS, n_cols=N_COLS):
    idx = np.arange(n_rows * n_cols).reshape(n_rows, n_cols)
    i = idx.take(range(idx.shape[axis] - k), axis=axis).ravel()
    j = idx.take(range(k, idx.shape[axis]), axis=axis).ravel()
    return list(zip(i, j))


def centroid(c, rad):
    """Intensity-weighted centroid of the surface within +-rad px of zero lag."""
    half = c.shape[0] // 2
    sl = slice(half - rad, half + rad + 1)
    w = np.clip(c[sl, sl], 0, None)
    y, x = np.mgrid[-rad:rad + 1, -rad:rad + 1]
    tot = w.sum()
    return float((w * x).sum() / tot), float((w * y).sum() / tot)


def surface_stats(c, ref_peak):
    """Peak height (relative to a k=0 reference) and centroids at two window radii.

    The two radii are a deliberate cross-check: a genuine translation moves the
    whole surface, so both centroids agree. A centroid driven by the broad,
    anisotropically-decaying tails collapses when the window is tightened.
    """
    half = c.shape[0] // 2
    cx4, cy4 = centroid(c, min(4, half))
    cx, cy = centroid(c, half)
    return dict(peak=c[half, half] / ref_peak, cx=cx, cy=cy, cx4=cx4, cy4=cy4)


def analyse(path, tag):
    with h5py.File(path, "r") as f:
        r = process(find_frames(f), N_ROWS, N_COLS, 64)
    F, shape, S = residual_ffts(r["crops"], r["roi_good"])
    surfaces, stats = {}, {}
    # k = 0: autocorrelation. Its zero lag is 1 + (shot noise)/(signal power), i.e. > 1 -- the
    # noise spike is the ONLY thing there that is not speckle, which is why it gets masked.
    surfaces[("auto", 0)] = mean_xcorr(
        F, shape, [(i, i) for i in range(len(F))], S, good=r["roi_good"])
    half = surfaces[("auto", 0)].shape[0] // 2
    for axis, name in ((1, "fast"), (0, "slow")):
        for k in (1, 2, 4):
            surfaces[(name, k)] = mean_xcorr(
                F, shape, neighbour_pairs(k, axis), S, good=r["roi_good"])
    for key in list(surfaces):
        if key[0] != "auto":
            stats[key] = surface_stats(surfaces[key], 1.0)
    print(f"{tag}:")
    for key in sorted(stats, key=lambda t: (t[0], t[1])):
        s = stats[key]
        print(f"   {key[0]:4s} k={key[1]} ({key[1] * STEP_NM:3.0f} nm): peak {s['peak']:.3f}  "
              f"centroid +-14px ({s['cx']:+.2f},{s['cy']:+.2f})  "
              f"+-4px ({s['cx4']:+.2f},{s['cy4']:+.2f})")
    return surfaces, stats


def main():
    fig, ax = plt.subplots(2, 4, figsize=(16, 10))
    for row, (path, tag, *_rest) in enumerate(SCANS):
        surfaces, stats = analyse(path, tag)
        auto = surfaces[("auto", 0)]
        half = auto.shape[0] // 2
        show = auto.copy()
        show[half, half] = np.nan                     # shot-noise delta lives only here
        panels = ((np.clip(show, 0, 1), "k=0 AUTOcorrelation\n(zero lag = shot noise, masked)"),
                  (surfaces[("fast", 1)], "CROSS, fast neighbour\nk=1 (20 nm)"),
                  (surfaces[("slow", 1)], "CROSS, slow neighbour\nk=1 (20 nm)"),
                  (surfaces[("slow", 4)], "CROSS, slow neighbour\nk=4 (80 nm)"))
        for col, (m, name) in enumerate(panels):
            im = ax[row, col].imshow(m, cmap="magma", vmin=0, vmax=1.0,
                                     extent=(-half, half, half, -half))
            ax[row, col].set_title(f"{tag}  {name}", fontsize=8.5)
            ax[row, col].set_xlabel("lag, det x [px]", fontsize=8)
            ax[row, col].set_ylabel("lag, det y [px]", fontsize=8)
            plt.colorbar(im, ax=ax[row, col], shrink=0.8).set_label("correlation coefficient", fontsize=7)
        for col, key in ((1, ("fast", 1)), (2, ("slow", 1)), (3, ("slow", 4))):
            s = stats[key]
            ax[row, col].text(0.03, 0.03, f"peak {s['peak']:.2f}\ncentroid "
                              f"({s['cx']:+.2f},{s['cy']:+.2f})", transform=ax[row, col].transAxes,
                              fontsize=7.5, color="w", va="bottom",
                              bbox=dict(fc="k", alpha=0.5, pad=2, ec="none"))
    fig.suptitle(
        "2D cross-correlation of neighbouring speckle IMAGES  C_k(dq) = <r_i(q) r_j(q+dq)>, "
        "envelope-subtracted\n"
        "k=0 is the autocorrelation (speckle size/shape, plus the shot-noise delta at zero lag); "
        "k>=1 are cross-correlations between scan neighbours.\n"
        "Peaks are TRUE correlation coefficients: each pair is normalised by sqrt(S_i S_j) with the "
        "Poisson term (S = sum r^2 - sum counts) removed, so k=0\n"
        "would be exactly 1 were its centre not corrupted by the shot-noise delta - the one thing "
        "there that is not speckle. The surfaces keep their shape and\n"
        "fall in height with k: the speckle DECORRELATES\n"
        "essentially in place. Centroids are quoted at two window radii as a check - a real "
        "translation is window-independent (the fast axis, ~0.08 px/step,\n"
        "linear in k) while a tail-driven artefact halves when the window tightens (the slow axis). "
        "Both are in tension with the COM field and are NOT yet physics.", fontsize=9)
    fig.tight_layout(rect=(0, 0, 1, 0.86), h_pad=3.0)
    fig.savefig("xcorr_2d.png", dpi=140)
    print("wrote xcorr_2d.png")


def selftest():
    """Known speckle + known decorrelation: check peak height, shape and centroid behave."""
    rng = np.random.default_rng(0)
    n, grain = 128, 3.0
    # Opposite-edge pixels are far apart, not one pixel apart: correlation must be linear,
    # without the periodic wraparound produced by an unpadded FFT.
    edge = np.zeros((2, n, n))
    edge[0, 0, 0], edge[1, -1, -1] = 1, 1
    Fe = np.array([np.fft.rfft2(f, s=(2 * n, 2 * n)) for f in edge])
    ce = mean_xcorr(Fe, (n, n), [(0, 1)], np.ones(2), half=2)
    assert np.max(np.abs(ce)) < 1e-12, "cross-correlation must not wrap across ROI edges"
    flat = np.ones((n, n))
    Fflat = np.array([np.fft.rfft2(flat, s=(2 * n, 2 * n))])
    cflat = mean_xcorr(Fflat, (n, n), [(0, 0)], np.array([flat.size]), half=2)
    assert np.allclose(cflat, 1, atol=1e-12), "finite overlap must not taper correlation"

    ky, kx = np.meshgrid(np.fft.fftfreq(n), np.fft.fftfreq(n), indexing="ij")
    filt = np.exp(-2 * (np.pi * grain) ** 2 * (ky ** 2 + kx ** 2) / 2)
    spk = lambda: np.fft.ifft2(np.fft.fft2(rng.standard_normal((n, n))) * filt).real
    a = spk()
    frames, mix = [], 0.6                       # 60% shared structure between "neighbours"
    for _ in range(40):
        frames.append(a)
        frames.append(mix * a + np.sqrt(1 - mix ** 2) * spk())
    F = np.array([np.fft.rfft2(f, s=(2 * n, 2 * n)).astype(np.complex64) for f in frames])
    pairs_same = [(2 * i, 2 * i) for i in range(40)]
    pairs_mix = [(2 * i, 2 * i + 1) for i in range(40)]
    S = np.array([(f ** 2).sum() for f in frames])
    c0 = mean_xcorr(F, (n, n), pairs_same, S)
    c1 = mean_xcorr(F, (n, n), pairs_mix, S)
    h = c0.shape[0] // 2
    got = c1[h, h] / c0[h, h]
    st = surface_stats(c1, c0[h, h])
    print(f"SELFTEST: imposed correlation {mix:.2f} -> recovered peak ratio {got:.3f}; "
          f"centroid ({st['cx']:+.3f}, {st['cy']:+.3f}) px")
    assert abs(got - mix) < 0.05, f"peak height should recover the imposed correlation: {got:.3f}"
    # 0.12 px: finite-sample scatter of the centroid with 40 pairs, still far below the 1 px
    # translation the next assertion has to detect
    assert abs(st["cx"]) < 0.12 and abs(st["cy"]) < 0.12, "unshifted speckle must stay centred"
    # a genuine translation must move the centroid, so the statistic can detect one
    shifted = [np.fft.ifft2(np.fft.fft2(f) * np.exp(-2j * np.pi * kx * 1.0)).real for f in frames]
    Fs = np.array([np.fft.rfft2(f, s=(2 * n, 2 * n)).astype(np.complex64) for f in shifted])
    Fm = np.array([F[i] if i % 2 == 0 else Fs[i] for i in range(len(F))])
    c2 = mean_xcorr(Fm, (n, n), pairs_mix, S)
    st2 = surface_stats(c2, c0[h, h])
    print(f"SELFTEST: imposed 1.0 px translation -> centroid ({st2['cx']:+.3f}, {st2['cy']:+.3f}) px")
    assert abs(st2["cx"] - 1.0) < 0.25, "an imposed translation must show in the centroid"
    print("SELFTEST PASS")


if __name__ == "__main__":
    p = argparse.ArgumentParser()
    p.add_argument("--selftest", action="store_true")
    a = p.parse_args()
    selftest() if a.selftest else main()


# --- the OTHER 2D slice: correlation vs SCAN offset, at zero detector lag -------------
# C(0; d_fast, d_slow). The panels above hold the scan offset fixed and show the detector
# lag plane; this holds the detector lag at zero and shows the SCAN plane, so the
# decorrelation anisotropy appears directly as an ellipse in real space.

def scan_space_corr(crops, good, max_off=6):
    """Poisson-corrected speckle correlation for every (d_fast, d_slow) scan offset.

    Done as one Gram matrix of normalised residuals: every frame pair at once, then
    averaged by scan offset. Cheap enough (2911 x 16k) and exact.
    """
    m = good.ravel()
    mean_pat = crops.mean(axis=0) * good
    flat_mean = mean_pat.ravel()[m]
    R = crops.reshape(len(crops), -1)[:, m].astype(np.float32)
    counts = R.sum(axis=1)
    R = R - np.outer(counts / flat_mean.sum(), flat_mean)
    S = np.maximum((R ** 2).sum(axis=1) - counts, 1e-9)
    R /= np.sqrt(S)[:, None].astype(np.float32)
    G = R @ R.T
    idx = np.arange(len(crops)).reshape(N_ROWS, N_COLS)
    n = 2 * max_off + 1
    out = np.full((n, n), np.nan)
    for si, ds in enumerate(range(-max_off, max_off + 1)):
        for fi, df in enumerate(range(-max_off, max_off + 1)):
            r0, r1 = max(0, -ds), min(N_ROWS, N_ROWS - ds)
            c0, c1 = max(0, -df), min(N_COLS, N_COLS - df)
            if r1 <= r0 or c1 <= c0:
                continue
            a = idx[r0:r1, c0:c1].ravel()
            b = idx[r0 + ds:r1 + ds, c0 + df:c1 + df].ravel()
            out[si, fi] = np.median(G[a, b])
    return out


def scan_space_figure():
    fig, ax = plt.subplots(1, 2, figsize=(11, 5))
    for col, (path, tag, *_rest) in enumerate(SCANS):
        with h5py.File(path, "r") as f:
            r = process(find_frames(f), N_ROWS, N_COLS, 64)
        C = scan_space_corr(r["crops"], r["roi_good"])
        h = C.shape[0] // 2
        ext = (-h * STEP_NM, h * STEP_NM, h * STEP_NM, -h * STEP_NM)
        im = ax[col].imshow(C, cmap="viridis", vmin=0, vmax=1, extent=ext)
        cs = ax[col].contour(np.linspace(ext[0], ext[1], C.shape[1]),
                             np.linspace(ext[3], ext[2], C.shape[0]), C,
                             levels=[0.3, 0.5, 0.7], colors="w", linewidths=0.8)
        ax[col].clabel(cs, fontsize=7, fmt="%.1f")
        ax[col].set_title(f"{tag}\nspeckle correlation vs SCAN offset", fontsize=10)
        ax[col].set_xlabel("fast (zpx) offset [nm]")
        ax[col].set_ylabel("slow (zpy) offset [nm]")
        plt.colorbar(im, ax=ax[col], shrink=0.85).set_label("correlation coefficient", fontsize=8)
        # half-correlation distance along each axis, by linear interpolation
        for axis, lab in ((1, "fast"), (0, "slow")):
            prof = C[h, h:] if axis == 1 else C[h:, h]
            d = np.arange(len(prof)) * STEP_NM
            below = np.where(prof < 0.5)[0]
            hw = np.interp(0.5, prof[below[0] - 1:below[0] + 1][::-1],
                           d[below[0] - 1:below[0] + 1][::-1]) if len(below) else np.nan
            print(f"{tag} {lab}: correlation falls to 0.5 at {hw:.0f} nm")
    fig.suptitle("The other 2D slice: C(0; d_fast, d_slow) - detector lag held at zero, SCAN plane shown.\n"
                 "Decorrelation anisotropy appears directly as an ellipse in real space "
                 "(elongated along fast = the 1/sin(theta)-stretched footprint direction).", fontsize=10)
    fig.tight_layout(rect=(0, 0.02, 1, 0.86))
    fig.savefig("xcorr_scanspace.png", dpi=140)
    print("wrote xcorr_scanspace.png")
