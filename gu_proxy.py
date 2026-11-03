#!/usr/bin/env python3
"""Rung 1: integrate the peak-COM field into a G.u phase-proxy surface.

The per-frame Bragg-peak COM shift measures the local mean q-offset, i.e. the
transverse gradient of the object phase phi = G.u (mixed lattice tilt + strain;
DPC / Shack-Hartmann logic). Least-squares integration (DCT Poisson solve,
Neumann BC) turns (COM_x, COM_y) into a smooth phase surface — the rung-1
deliverable and the initialization object for the rung-2 BPP solver, in the
spirit of Hruszkewycz PRL 110, 177601 (2013) picometre displacement maps.

Detector-to-sample axis pairing/signs are not yet fixed by a geometry layer, so
all 8 pairings are tried and the one minimizing the non-integrable residual
|g - grad(Phi)| / |g| is used (reported; the residual doubles as an
integrability diagnostic). Scale assumes Eiger2 75 um px at R = 522.5 mm,
E = 10.0 keV, 20 nm steps — pending David's confirmation, so pm values are
provisional.
"""
import numpy as np
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy.fft import dctn, idctn
from scipy.ndimage import gaussian_filter

from speckle_geometry import DIST_M, LAMBDA_NM, PX_M, STEP_NM

PX_RAD = PX_M / DIST_M           # detector angular pixel [rad] (DIST_M pending David: 4.5% on pm scale)
K = 2 * np.pi / (10 * LAMBDA_NM)  # [1/A]
STEP_A = 10 * STEP_NM            # sample-frame step [A]
RAD_PER_PX_STEP = K * PX_RAD * STEP_A  # phase gradient [rad/step] per COM px
G_HKL = {"202": 2 * np.pi / 1.55210, "131": 2 * np.pi / 1.32365}  # [1/A], a=4.39


def integrate(gx, gy):
    """Least-squares integrate gradient field -> (Phi, residual fraction)."""
    m, n = gx.shape
    div = np.zeros((m, n))
    div[:, 1:] += gx[:, 1:] - gx[:, :-1]
    div[1:, :] += gy[1:, :] - gy[:-1, :]
    lam = (2 * np.cos(np.pi * np.arange(m) / m) - 2)[:, None] + \
          (2 * np.cos(np.pi * np.arange(n) / n) - 2)[None, :]
    lam[0, 0] = 1.0
    phat = dctn(div, norm="ortho") / lam
    phat[0, 0] = 0.0
    phi = idctn(phat, norm="ortho")
    rx, ry = np.diff(phi, axis=1) - gx[:, :-1], np.diff(phi, axis=0) - gy[:-1, :]
    res = np.sqrt(np.mean(rx ** 2) + np.mean(ry ** 2)) / \
          np.sqrt(np.mean(gx ** 2) + np.mean(gy ** 2))
    return phi, float(res)


def decompose(comx, comy, sigma=2.0):
    """Helmholtz-decompose the low-passed COM field at the most-gradient frame.

    The detector->sample frame angle is unknown and mixes div/curl, so scan a
    pure rotation (+reflection) and keep the angle maximizing the
    gradient-compatible power: grad_frac = 1 - res^2 is then an UPPER BOUND on
    how much of the domain-scale COM field can be any scalar strain-gradient
    G.u surface; the remainder is irreducibly rotational (twist). At that
    frame, Phi (from div, DCT/Neumann) is the envelope G.u proxy and Psi (from
    curl, same solver — Neumann is approximate for the rotational part, fine
    for a display map) is the twist stream-function whose domains are the
    real-space twist texture. sigma low-pass isolates the domain envelope.
    """
    # mean (not median) subtraction: a constant COM offset is frame calibration,
    # and the Neumann solver cannot represent it, so it must not enter the field
    a = (comx - np.mean(comx)) * RAD_PER_PX_STEP
    b = (comy - np.mean(comy)) * RAD_PER_PX_STEP
    a_s = gaussian_filter(a, sigma) - np.mean(a)
    b_s = gaussian_filter(b, sigma) - np.mean(b)
    best = None
    for refl in (1, -1):
        for al in np.deg2rad(np.arange(0, 360, 2)):
            gx = np.cos(al) * a_s - np.sin(al) * refl * b_s
            gy = np.sin(al) * a_s + np.cos(al) * refl * b_s
            phi, res = integrate(gx, gy)
            if best is None or res < best[1]:
                best = (phi, res, (np.rad2deg(al), refl), (gx, gy))
    phi, res, (alpha, refl), (gx, gy) = best
    # twist stream-function from the curl of the same field (rotated 90deg trick:
    # curl(gx,gy) = div(gy,-gx), so reuse the Poisson solver); remove its LSQ
    # plane — a planar Psi is a constant rotational field, i.e. residual frame
    # offset, not twist texture
    psi, _ = integrate(gy, -gx)
    yy, xx = np.mgrid[0:psi.shape[0], 0:psi.shape[1]]
    A = np.stack([np.ones(psi.size), xx.ravel(), yy.ravel()], axis=1)
    psi = psi - (A @ np.linalg.lstsq(A, psi.ravel(), rcond=None)[0]).reshape(psi.shape)
    return phi, psi, float(1 - res ** 2), (alpha, refl)


def main():
    fig, ax = plt.subplots(2, 4, figsize=(16, 9))
    fracs = {}
    for row, (scan, hkl) in enumerate((("data/scan_309_100309_maps.npz", "202"),
                                       ("data/scan_313_100313_maps.npz", "131"))):
        Z = np.load(scan)
        phi, psi, grad_frac, (alpha, refl) = decompose(Z["comx"], Z["comy"])
        fracs[hkl] = grad_frac
        tw_int = np.corrcoef(psi.ravel(), Z["total"].ravel())[0, 1]
        u_pm = phi / G_HKL[hkl] * 100.0  # A -> pm
        tag = f"({hkl})"
        for col, (m, name, cmap, unit) in enumerate((
                (Z["comx"] - np.median(Z["comx"]), "COM_x (input)", "RdBu_r", "det px"),
                (Z["comy"] - np.median(Z["comy"]), "COM_y (input)", "RdBu_r", "det px"),
                (u_pm, f"gradient part: G·u proxy\n(<= {100 * grad_frac:.0f}% of field power)",
                 "RdBu_r", "u_G [pm] (provisional scale)"),
                (psi, f"rotational part: twist map Psi\n(>= {100 * (1 - grad_frac):.0f}% of field power)",
                 "PuOr_r", "stream fn [arb]"))):
            v = np.nanpercentile(np.abs(m), 98)
            im = ax[row, col].imshow(m, cmap=cmap, vmin=-v, vmax=v)
            ax[row, col].set_title(f"{tag} {name}", fontsize=9)
            ax[row, col].set_xlabel("fast step", fontsize=8)
            ax[row, col].set_ylabel("slow step", fontsize=8)
            plt.colorbar(im, ax=ax[row, col], shrink=0.8).set_label(unit, fontsize=7)
        print(f"{tag}: best frame alpha={alpha:.0f}deg refl={refl}, gradient-compatible <= {grad_frac:.2f}, "
              f"u_G span {np.nanpercentile(u_pm, 1):.0f}..{np.nanpercentile(u_pm, 99):.0f} pm, "
              f"corr(twist, intensity) = {tw_int:+.2f}")
        assert 0.0 < grad_frac < 1.0 and np.isfinite(u_pm).all()
    fig.suptitle(
        "Rung 1: Helmholtz decomposition of the domain-scale peak-COM field (2-step low-pass), "
        "at the most-gradient-compatible detector->sample frame\n"
        f"RESULT: <={100 * fracs['202']:.0f}%/{100 * fracs['131']:.0f}% of the (202)/(131) COM field can be "
        "a scalar strain-gradient G·u surface - the majority is rotational: real-space TWIST-dominated mosaic\n"
        "(nanoscale counterpart of the ensemble twist mosaicity; quantitative case for 3-angle tilt/strain "
        "separation + phasing over COM shortcuts). pm scale assumes Eiger2 75 um, 10 keV, 20 nm steps, "
        "R = 522.5 mm readback (= 500 mm nominal setting: 4.5% scale ambiguity, pending David)", fontsize=9)
    fig.tight_layout(rect=(0, 0, 1, 0.92))
    fig.savefig("gu_proxy.png", dpi=140)
    print("wrote gu_proxy.png")


if __name__ == "__main__":
    main()
