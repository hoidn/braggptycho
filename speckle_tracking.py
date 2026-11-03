#!/usr/bin/env python3
"""Bragg speckle tracking: 2D cross-correlation of adjacent patterns -- and why it fails here.

This implements the senior PI's meeting-1 proposal in full. The zero-lag value
of the cross-correlation is the similarity number quicklook already maps; the
PEAK POSITION is the extra information -- how far the speckle translated on the
detector when the probe stepped, i.e. what David watched happen live. Two
neighbours (right, below) x two components = literally his "4 numbers for each
image".

RESULT: on the shared scans it does NOT work, and the reason is diagnostic
rather than incidental -- 75-82% of the fitted sub-pixel shifts land within
0.05 px of an integer lag, where an unbiased estimator gives 10%. That is
PEAK LOCKING: the cross-correlation peak is only ~1.3 px wide, so parabolic
sub-pixel refinement is pinned to integer lags. The exception proves it --
the one component aligned with the BROAD speckle axis is unbiased (14.6%).

So the failure is the direct consequence of the speckle undersampling measured
in speckle_geometry.py (1.3-1.7 px FWHM, below the 2 px Nyquist floor), and it
converts that abstract sampling number into a concrete cost: an entire analysis
modality is unavailable at this detector distance. Move the detector back (or
use more symmetric reflections) and this becomes a working measurement -- it is
the strongest argument for the acquisition change, and a next-beamtime ask.

Corroborating checks, all negative in a consistent way:
  - shift increments do not correlate with COM increments (|r| < 0.08, all
    components, both scans);
  - larger strides make it worse, not better: by offset 8 the "shift" rms
    reaches 3.1 px against a true COM change of 0.5 px, i.e. the argmax is
    jumping randomly inside the search window once the speckle decorrelates
    (envelope-subtracted speckle decorrelates over 82-161 nm, so 160 nm apart
    there is nothing left to track);
  - the measured shift rms (0.031 px) is SMALLER than the real signal
    (~0.08-0.12 px from COM), i.e. systematic under-response, not scatter.

Precedent: X-ray speckle tracking (Morgan JAC 2020) is established in
TRANSMISSION; no Bragg implementation exists. The gap is real, but this dataset
cannot fill it.

Usage: python speckle_tracking.py [--selftest]
"""
import argparse

import numpy as np
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt

from quicklook import find_frames, process
from speckle_geometry import SCANS

try:
    import hdf5plugin  # noqa: F401  (registers bslz4)
except ImportError:
    pass
import h5py

N_ROWS, N_COLS = 71, 41
LOCK_TOL = 0.05          # |shift - nearest integer| counted as "locked"
LOCK_UNIFORM = 2 * LOCK_TOL  # unbiased expectation for that fraction


def residuals(crops, good):
    """Envelope-subtracted frames: the speckle alone, flux-matched to the scan mean."""
    mean_pat = crops.mean(axis=0) * good
    denom = mean_pat.sum()
    out = np.empty_like(crops)
    for i, f in enumerate(crops):
        r = f * good
        out[i] = r - (r.sum() / denom) * mean_pat
    return out


def xcorr_shift(a, b, max_lag=6):
    """Sub-pixel shift of `b` relative to `a` (dy, dx) from the 2D cross-correlation peak.

    Parabolic refinement about the integer peak -- the standard estimator, and
    the one that peak-locks when the correlation peak is ~1 px wide.
    """
    A, B = np.fft.rfft2(a), np.fft.rfft2(b)
    c = np.fft.fftshift(np.fft.irfft2(np.conj(A) * B, s=a.shape))
    cy, cx = np.array(a.shape) // 2
    win = c[cy - max_lag:cy + max_lag + 1, cx - max_lag:cx + max_lag + 1]
    iy, ix = np.unravel_index(np.argmax(win), win.shape)
    iy, ix = iy + cy - max_lag, ix + cx - max_lag
    if not (0 < iy < c.shape[0] - 1 and 0 < ix < c.shape[1] - 1):
        return np.nan, np.nan, c[iy, ix]

    def refine(m, p):
        d = m[0] - 2 * m[1] + m[2]
        return p + (0.5 * (m[0] - m[2]) / d if abs(d) > 1e-12 else 0.0)

    dy = refine(c[iy - 1:iy + 2, ix], iy) - cy
    dx = refine(c[iy, ix - 1:ix + 2], ix) - cx
    return dy, dx, c[iy, ix]


def track(res, n_rows=N_ROWS, n_cols=N_COLS, max_lag=6):
    """Speckle shift to the next frame along each scan axis -- the '4 numbers per image'."""
    idx = np.arange(len(res)).reshape(n_rows, n_cols)
    out = {k: np.full((n_rows, n_cols), np.nan) for k in
           ("fast_dy", "fast_dx", "slow_dy", "slow_dx")}
    for r in range(n_rows):
        for c in range(n_cols):
            if c + 1 < n_cols:
                dy, dx, _ = xcorr_shift(res[idx[r, c]], res[idx[r, c + 1]], max_lag)
                out["fast_dy"][r, c], out["fast_dx"][r, c] = dy, dx
            if r + 1 < n_rows:
                dy, dx, _ = xcorr_shift(res[idx[r, c]], res[idx[r + 1, c]], max_lag)
                out["slow_dy"][r, c], out["slow_dx"][r, c] = dy, dx
    return out


def lock_fraction(v):
    """Fraction of sub-pixel estimates sitting within LOCK_TOL of an integer lag.

    LOCK_UNIFORM (10%) for an unbiased estimator; >> that means peak locking.
    """
    v = np.asarray(v)
    v = v[np.isfinite(v)]
    return float(np.mean(np.abs(v - np.round(v)) < LOCK_TOL))


def analyse(path, tag):
    with h5py.File(path, "r") as f:
        r = process(find_frames(f), N_ROWS, N_COLS, 64)
    sh = track(residuals(r["crops"], r["roi_good"]))
    locks = {k: lock_fraction(v) for k, v in sh.items()}
    # increment-vs-increment consistency (no integration, so no random-walk build-up)
    dcom = {"fast_dx": np.diff(r["comx"], axis=1), "fast_dy": np.diff(r["comy"], axis=1),
            "slow_dx": np.diff(r["comx"], axis=0), "slow_dy": np.diff(r["comy"], axis=0)}
    print(f"{tag}:")
    for k in ("fast_dx", "fast_dy", "slow_dx", "slow_dy"):
        s = sh[k][:, :-1] if k.startswith("fast") else sh[k][:-1, :]
        a, b = s.ravel(), dcom[k].ravel()
        m = np.isfinite(a) & np.isfinite(b)
        print(f"   {k}: locked {100 * locks[k]:5.1f}% (unbiased {100 * LOCK_UNIFORM:.0f}%)  "
              f"rms {np.nanstd(a):.3f} px vs dCOM {np.nanstd(b):.3f} px  "
              f"corr {np.corrcoef(a[m], b[m])[0, 1]:+.3f}")
    return sh, locks


def main():
    fig, ax = plt.subplots(2, 3, figsize=(14, 8))
    for row, (path, tag, *_rest) in enumerate(SCANS):
        sh, locks = analyse(path, tag)
        best = min(locks, key=lambda k: locks[k])      # least peak-locked component
        worst = max(locks, key=lambda k: locks[k])
        for col, key, lab in ((0, best, "least peak-locked"), (1, worst, "most peak-locked")):
            m = sh[key] - np.nanmedian(sh[key])
            v = np.nanpercentile(np.abs(m[np.isfinite(m)]), 96)
            im = ax[row, col].imshow(m, cmap="RdBu_r", vmin=-v, vmax=v)
            ax[row, col].set_title(f"{tag}  {key}\n{lab}: {100 * locks[key]:.0f}% at integer lags",
                                   fontsize=9)
            ax[row, col].set_xlabel("fast step", fontsize=8)
            ax[row, col].set_ylabel("slow step", fontsize=8)
            plt.colorbar(im, ax=ax[row, col], shrink=0.8).set_label("det px / step", fontsize=7)
        for k, c in zip(("fast_dx", "fast_dy", "slow_dx", "slow_dy"), ("C0", "C2", "C1", "C3")):
            v = sh[k][np.isfinite(sh[k])]
            ax[row, 2].hist(np.abs(v - np.round(v)), bins=20, range=(0, 0.5), histtype="step",
                            density=True, color=c, lw=1.4,
                            label=f"{k}: {100 * locks[k]:.0f}% locked")
        ax[row, 2].axhline(1.0, color="k", ls="--", lw=1, label="unbiased (flat)")
        ax[row, 2].set_title(f"{tag}  distance to nearest integer lag", fontsize=9)
        ax[row, 2].set_xlabel("|shift - nearest integer| [px]", fontsize=8)
        ax[row, 2].set_ylabel("density", fontsize=8)
        ax[row, 2].legend(fontsize=6.5)
    fig.suptitle(
        "Bragg speckle tracking (the PI's \"4 numbers per image\"): implemented, and diagnosed as "
        "BLOCKED by detector sampling\n"
        "The sub-pixel cross-correlation shifts pile up at integer lags (75-82% within 0.05 px, vs 10% "
        "unbiased) = PEAK LOCKING, because the\n"
        "correlation peak is only ~1.3 px wide. The one component aligned with the BROAD speckle axis "
        "is unbiased (15%), which proves the mechanism.\n"
        "=> This is the concrete cost of the undersampling measured in speckle_geometry.py: a whole "
        "analysis modality is unavailable at this detector\n"
        "distance. Estimator itself is sound (selftest recovers imposed shifts to 0.011 px). Move the "
        "detector back -> the measurement works.", fontsize=9)
    fig.tight_layout(rect=(0, 0, 1, 0.87))
    fig.savefig("speckle_tracking.png", dpi=140)
    print("wrote speckle_tracking.png")


def _speckle(rng, n, grain_px):
    """Synthetic speckle field whose correlation length is ~grain_px pixels."""
    ky, kx = np.meshgrid(np.fft.fftfreq(n), np.fft.fftfreq(n), indexing="ij")
    f = np.exp(-2 * (np.pi * grain_px) ** 2 * (ky ** 2 + kx ** 2) / 2)
    return np.fft.ifft2(np.fft.fft2(rng.standard_normal((n, n))) * f).real


def selftest():
    """(1) the estimator is accurate on well-sampled speckle; (2) it peak-locks when undersampled."""
    rng = np.random.default_rng(0)
    n = 128
    ky, kx = np.meshgrid(np.fft.fftfreq(n), np.fft.fftfreq(n), indexing="ij")
    base = _speckle(rng, n, 3.0)
    truth = [(0.0, 0.0), (0.37, -0.62), (1.5, 2.0), (-2.25, 0.85)]
    errs = []
    for ty, tx in truth:
        shifted = np.fft.ifft2(np.fft.fft2(base) * np.exp(-2j * np.pi * (ky * ty + kx * tx))).real
        dy, dx, _ = xcorr_shift(base, shifted)
        errs.append(max(abs(dy - ty), abs(dx - tx)))
    worst = max(errs)
    print(f"SELFTEST 1 (well-sampled speckle): worst sub-pixel error {worst:.3f} px")
    assert worst < 0.08, f"sub-pixel shift recovery off by {worst:.3f} px"
    dy, dx, _ = xcorr_shift(base, base * 1.7)
    assert abs(dy) < 1e-6 and abs(dx) < 1e-6, "intensity scaling leaked into the shift"

    # SELFTEST 2: reproduce the failure mechanism seen on the real data. Random sub-pixel
    # shifts applied to WELL-sampled vs UNDERSAMPLED speckle; only the latter must peak-lock.
    locks = {}
    for lab, grain in (("well-sampled", 3.0), ("undersampled", 0.55)):
        got = []
        for _ in range(200):
            fld = _speckle(rng, n, grain)
            ty, tx = rng.uniform(-0.5, 0.5, 2)
            sh = np.fft.ifft2(np.fft.fft2(fld) * np.exp(-2j * np.pi * (ky * ty + kx * tx))).real
            d = xcorr_shift(fld, sh)
            got.extend([d[0], d[1]])
        locks[lab] = lock_fraction(got)
        print(f"SELFTEST 2 ({lab}, grain {grain} px): {100 * locks[lab]:.0f}% of estimates at "
              f"integer lags (unbiased {100 * LOCK_UNIFORM:.0f}%)")
    assert locks["well-sampled"] < 0.25, "well-sampled speckle must not peak-lock"
    assert locks["undersampled"] > 2 * locks["well-sampled"], \
        "undersampled speckle must peak-lock -- this is the diagnosis made on the real data"
    print("SELFTEST PASS")


if __name__ == "__main__":
    p = argparse.ArgumentParser()
    p.add_argument("--selftest", action="store_true")
    a = p.parse_args()
    selftest() if a.selftest else main()
