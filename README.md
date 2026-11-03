# braggptycho

Analysis of Bragg nanodiffraction scans of epitaxial δ-NbN(111) thin films on sapphire, from the 26-ID nanoprobe. Nothing here reconstructs an image. The scripts compute maps and statistics straight from the stacks of diffraction frames: intensity, Bragg-peak centre of mass, speckle correlation, speckle size, and photon budget.

There are two real scans, both from a 10 nm film:

| Scan | Reflection | θ | Raster |
|---|---|---|---|
| `scan_309_100309` | NbN (202) | 11.2° | 71 × 41 (slow zpy × fast zpx), 20 nm steps |
| `scan_313_100313` | NbN (131) | 18.6° | 71 × 41 |

`BRIEFING.md` has the science context, the verified geometry (E ≈ 10 keV, Eiger2 75 µm pixels, R = 522.5 mm) and the findings. `SESSION_LOG_2026-08-20.md` records how the analysis was done.

## Key figures

**Speckle cross-correlation vs detector lag** (`xcorr_2d.py`). k = 0 is the autocorrelation. For k ≥ 1, each surface correlates a frame with the frame k scan steps away. The peak height gives the decorrelation, the width and shape give the coherently lit extent, and the centroid gives any systematic speckle translation.

![xcorr_2d](xcorr_2d.png)

**Speckle correlation vs scan offset, at zero detector lag** (`xcorr_2d.scan_space_figure()`). The decorrelation anisotropy appears as an ellipse in real space.

![xcorr_scanspace](xcorr_scanspace.png)

Other figures: [twin comparison](twin_comparison.png) · [speckle geometry](speckle_geometry.png) · [speckle tracking](speckle_tracking.png) · [G·u phase proxy](gu_proxy.png) · [photon budget](photon_budget.png) · [rocking-curve Fisher information](fisher_rocking.png) · quicklooks for [scan309](data/scan_309_100309_quicklook.png) and [scan313](data/scan_313_100313_quicklook.png) · [synthetic selftest scan](synthetic_scan_quicklook.png)

## Data

- **Raw `.h5` files (about 700 MB each) are not in the repo.** Each file holds one detector frame per scan point in `/entry/data/data`, compressed with bitshuffle/LZ4, so `hdf5plugin` is required.
- **`data/*_maps.npz`** are the reduced maps that `quicklook.py` wrote with default flags. The scripts that read only these maps run without the raw files.
- **`synthetic_scan.h5`** is the synthetic twin-domain scan written by `python quicklook.py --selftest`.

## Scripts

| Script | Input | Output | What it shows |
|---|---|---|---|
| `quicklook.py <scan>.h5` | raw h5 | `<scan>_maps.npz`, `<scan>_quicklook.png` (saved next to the input) | Masks bad pixels and infers the raster grid. Plots the mean pattern, total and max-pixel maps, peak COM, adjacent-frame correlation, and correlation decay. |
| `compare_twins.py` | `data/*_maps.npz` | `twin_comparison.png` | (202) vs (131) side by side. The raw correlation differs between them; after Poisson correction they agree. |
| `gu_proxy.py` | `data/*_maps.npz` | `gu_proxy.png` | Integrates the peak-COM field into a G·u phase-proxy surface. The pm scale is provisional. |
| `photon_budget.py` | `data/*_maps.npz` | `photon_budget.png` | Counts per pixel vs radius from the Bragg peak. |
| `fisher_rocking.py` | none (analytic) | `fisher_rocking.png` | Where on the rocking curve one photon carries the most strain information. |
| `speckle_geometry.py` | raw h5 | `speckle_geometry.png` | Fits the speckle autocorrelation to get detector sampling and the extents of the coherent illumination. Defines `SCANS`, which the other raw-h5 scripts import. |
| `speckle_tracking.py` | raw h5 | `speckle_tracking.png` | Tracks how far the speckle shifts between neighbouring frames. It fails on this data because of peak locking, a consequence of speckle undersampling. |
| `xcorr_2d.py` | raw h5 | `xcorr_2d.png` | 2D speckle cross-correlation surfaces vs detector lag, for scan neighbours k = 0 (autocorrelation) through k ≥ 1. |
| `xcorr_2d.scan_space_figure()` | raw h5 | `xcorr_scanspace.png` | Speckle correlation vs scan offset at zero detector lag. Not run by the CLI; call it from Python (see below). |

Every script that has `--selftest` validates itself on synthetic data and needs no real files: `quicklook.py`, `speckle_geometry.py`, `speckle_tracking.py`, `xcorr_2d.py`. The first, second and fourth pass. `speckle_tracking.py` fails selftest 2: undersampled synthetic speckle (grain 0.55 px) puts only 12% of shift estimates at integer lags, while the assertion requires more than twice the well-sampled rate of 11%. So the synthetic model does not yet reproduce the peak locking the script reports on the real data.

## Running

The scripts need `numpy`, `scipy`, `matplotlib`, `h5py` and `hdf5plugin`. Run them from the repo root, because input paths are relative.

```bash
python quicklook.py --selftest              # self-checks, no data needed
python compare_twins.py                     # uses data/*_maps.npz

# scripts that read the raw h5: first set SCRATCH in speckle_geometry.py
# to the directory that holds scan_309_100309.h5 and scan_313_100313.h5
python xcorr_2d.py
python -c "import xcorr_2d; xcorr_2d.scan_space_figure()"
```

**Known issue:** `SCRATCH` in `speckle_geometry.py` still points to a temporary directory from the original session, and that directory no longer exists. Change it before running any script that reads the raw h5.
