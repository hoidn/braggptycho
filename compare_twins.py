#!/usr/bin/env python3
"""Meeting figure: (202) vs (131) twin reflections, reconstruction-free maps.

Reads the *_maps.npz written by quicklook.py (default flags) for scans 309
(202) and 313 (131) and renders a 2x4 comparison: integrated intensity,
peak COM, raw adjacent-frame speckle correlation, and Poisson-noise-corrected
correlation. The right-hand pair is the headline: raw correlation differs
(0.89 vs 0.65) but the corrected maps agree (~0.98-0.99), so the contrast is
photon statistics, not structure.
"""
import sys

import numpy as np
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt


def annotate(ax, m):
    ax.text(0.03, 0.03, f"median {np.nanmedian(m):.3f}", transform=ax.transAxes,
            fontsize=9, color="w", va="bottom",
            bbox=dict(fc="k", alpha=0.5, pad=2, ec="none"))


def main(a="data/scan_309_100309_maps.npz", b="data/scan_313_100313_maps.npz"):
    A, B = np.load(a), np.load(b)
    fig, ax = plt.subplots(2, 4, figsize=(16, 9))
    for row, (Z, tag) in enumerate([(A, "(202) twin A scan309"), (B, "(131) twin B scan313")]):
        tot = Z["total"]
        im = ax[row, 0].imshow(tot, cmap="viridis")
        ax[row, 0].set_title(f"{tag}\nintegrated intensity  (std/med {100 * tot.std() / np.median(tot):.1f}%)",
                             fontsize=10)
        plt.colorbar(im, ax=ax[row, 0], shrink=0.8).set_label("counts/frame", fontsize=8)
        # centre COM colormap on its median for tilt/strain contrast
        c = Z["comy"] - np.nanmedian(Z["comy"])
        v = np.nanpercentile(np.abs(c), 98)
        im = ax[row, 1].imshow(c, cmap="RdBu_r", vmin=-v, vmax=v)
        ax[row, 1].set_title(f"{tag}\npeak COM_y  (mixed tilt+strain)", fontsize=10)
        plt.colorbar(im, ax=ax[row, 1], shrink=0.8).set_label("det. pixel", fontsize=8)
        for col, key, name in ((2, "corr_y", "speckle corr (slow, raw)"),
                               (3, "corrc_y", "corr (slow, shot-noise corrected)")):
            im = ax[row, col].imshow(Z[key], cmap="inferno", vmin=0.5, vmax=1.0)
            ax[row, col].set_title(f"{tag}\n{name}", fontsize=10)
            annotate(ax[row, col], Z[key])
            plt.colorbar(im, ax=ax[row, col], shrink=0.8).set_label("Pearson r", fontsize=8)
        for a_ in ax[row]:
            a_.set_xlabel("fast (zpx) step", fontsize=8)
            a_.set_ylabel("slow (zpy) step", fontsize=8)
    fig.suptitle(
        "NbN 10 nm film - reconstruction-free maps, two twin reflections (DIFFERENT regions)\n"
        "scan309 (202) th=11.2deg / scan313 (131) th=18.6deg, E=10.0 keV (indexed from logged angles); "
        "grid 71x41, 20 nm sample-frame steps (sin-th prescaled zpx, beamline convention)\n"
        "resolution is probe-limited (~150 nm spot) and 1/sin(th)-elongated along x; "
        "right pair: identical corrected corr => raw contrast is pure photon statistics",
        fontsize=11)
    fig.tight_layout(rect=(0, 0, 1, 0.92))
    out = "twin_comparison.png"
    fig.savefig(out, dpi=140)
    print("wrote", out)
    for Z, tag in [(A, "202"), (B, "131")]:
        tot = Z["total"]
        print(f"{tag}: I med={np.median(tot):.0f} std/med={100 * tot.std() / np.median(tot):.1f}% "
              f"raw corr x/y={Z['decay_x'][0]:.3f}/{Z['decay_y'][0]:.3f} "
              f"corrected={Z['decayc_x'][0]:.3f}/{Z['decayc_y'][0]:.3f}")


if __name__ == "__main__":
    main(*(sys.argv[1:3] if len(sys.argv) >= 3 else []))
