#!/usr/bin/env python3
"""Plan-B seed: the photon budget of the real scans (radial counts/pixel).

Reads the *_maps.npz (mean_pattern + roi) and plots median counts/pixel/frame
vs radius from the Bragg peak: the speckle tail runs at ones-to-tens of
photons, which is the operating regime any reconstruction (and PtychoPINN's
Poisson loss) must be designed for. Also the count target for forward sims.
"""
import numpy as np
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt

fig, ax = plt.subplots(figsize=(8, 5))
for path, tag, color in (("data/scan_309_100309_maps.npz", "(202) scan309", "C0"),
                         ("data/scan_313_100313_maps.npz", "(131) scan313", "C1")):
    Z = np.load(path)
    mp = Z["mean_pattern"]
    y0, y1, x0, x1 = Z["roi"]
    cy, cx = (y0 + y1) // 2, (x0 + x1) // 2
    yy, xx = np.mgrid[0:mp.shape[0], 0:mp.shape[1]]
    rr = np.hypot(yy - cy, xx - cx)
    sel = (rr < 400) & (mp > 0)  # mp==0 marks masked bad pixels / empty modules
    rbin = rr[sel].astype(int)
    prof = np.array([np.median(mp[sel][rbin == k]) for k in range(1, 400)
                     if np.any(rbin == k)])
    ax.semilogy(prof, color, lw=1.5, label=f"{tag}  (ROI median {np.median(mp[y0:y1, x0:x1]):.1f} ph/px)")
    assert np.median(prof[:30]) > np.median(prof[-60:]), "profile should fall off the peak"
ax.axhline(10, color="gray", ls="--", lw=1)
ax.axhline(1, color="gray", ls=":", lw=1)
ax.text(360, 12, "10 ph/px", fontsize=9, color="gray", ha="right")
ax.text(360, 1.2, "1 ph/px", fontsize=9, color="gray", ha="right")
ax.set_xlabel("radius from Bragg peak [detector px]")
ax.set_ylabel("median counts / pixel / frame  (0.5 s)")
ax.set_title("Photon budget of the shared 10 nm scans: speckle tails at ones-to-tens of photons\n"
             "=> Poisson-aware / prior-regularized inversion; sets counts for forward simulations",
             fontsize=11)
ax.legend(fontsize=9)
fig.tight_layout()
fig.savefig("figures/photon_budget.png", dpi=140)
print("wrote figures/photon_budget.png")
