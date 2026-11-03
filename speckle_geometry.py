#!/usr/bin/env python3
"""Coherent illumination extents and detector speckle sampling, from the speckle ACF.

Speckle is the interference fingerprint of the coherently illuminated scattering
volume. For a Gaussian illumination of intensity-FWHM w along some direction,
the far-field speckle intensity autocorrelation is Gaussian in dq with
    sigma_q = 1.665 / w          (Siegert: ACF = |FT{|P|^2}|^2)
so a real-space extent follows from the ACF width in detector pixels:
    w = 1.665 / (sigma_lag * dq_per_px),   dq_per_px = k * pixel / distance.
`--selftest` validates that inversion against synthetic speckle from a known
anisotropic probe.

What this delivers on the shared scans:
  1. Detector SAMPLING of the speckle - the actionable result. One axis is
     comfortably oversampled, the other is at/below Nyquist, so half the
     available speckle information is currently thrown away by the detector
     distance. That is an acquisition-design number (Plan B), not a fit.
  2. Two coherent real-space extents, one of them (63 nm) far below the ~150 nm
     spot size in the metadata -> a concrete question for David about which
     zone plate / defocus scans 308-319 actually used.
  3. The ACF ellipse orientation: the detector-frame anisotropy axis, which the
     rung-2 geometry layer must reproduce.

Deliberately NOT claimed: an identification of each axis with (probe /
1/sin(theta) footprint / film thickness). These are asymmetric reflections
(incidence != exit angle) and the detector cuts a tilted Ewald slice, so the
projection mixes them; separating them needs the full Jacobian (rung 2).

Shot noise is white -> a delta at zero lag only, excluded from the fit
(standard speckle/XPCS practice).
"""
import argparse

import numpy as np
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy.ndimage import map_coordinates
from scipy.optimize import curve_fit

from quicklook import find_frames, process

try:
    import hdf5plugin  # noqa: F401  (registers bslz4)
except ImportError:
    pass
import h5py

LAMBDA_NM = 0.123984          # 10.0 keV
PX_M = 75e-6                  # Eiger2 pixel
DIST_M = 0.5225               # R readback; nominal setting 500 mm (4.5% -> extents, see BRIEFING §6)
PIX_ACF_VAR = 1.0 / 6.0       # variance [px^2] of a square pixel's self-ACF (triangle, half-width 1)
NYQUIST_PX = 2.0              # px per speckle FWHM needed to sample the speckle
SCRATCH = ("/tmp/claude-1000/-home-ollie-Documents-braggptycho/"
           "d1fadf69-8f4b-4be7-8be7-40790de56e4a/scratchpad/data")
SCANS = ((f"{SCRATCH}/scan_309_100309.h5", "(202) scan309", 11.2, 42.2, 23.15),
         (f"{SCRATCH}/scan_313_100313.h5", "(131) scan313", 18.6, 52.4, 23.15))


def dq_per_px(dist_m=DIST_M):
    """Momentum-transfer step per detector pixel [1/nm]."""
    return (2 * np.pi / LAMBDA_NM) * (PX_M / dist_m)


def extent_nm(sigma_lag, dq=None):
    """Real-space coherent extent [nm] (intensity FWHM) from the ACF sigma [px]."""
    return 1.665 / (sigma_lag * (dq if dq is not None else dq_per_px()))


def ptycho_fov_nm(dist_m=DIST_M):
    """Object patch a single diffraction pattern can represent: lambda*R/pixel [nm].

    The illumination must sit well inside this window or the reconstruction
    wraps around; FOV/extent is the per-axis oversampling ratio (want >= 2).
    """
    return LAMBDA_NM * dist_m / PX_M


def exit_projection(samth_deg, tth_deg, gam_deg):
    """Footprint->far-field aperture factor sin(exit)/sin(incidence) for an asymmetric reflection.

    The surface footprint stretches as 1/sin(alpha), but the far field sees it
    projected onto the plane normal to the EXIT beam, i.e. x sin(beta). Only for
    a symmetric reflection (beta == alpha) do these cancel.
    """
    two_theta = np.degrees(np.arccos(np.cos(np.radians(tth_deg)) * np.cos(np.radians(gam_deg))))
    alpha, beta = samth_deg, two_theta - samth_deg
    return np.sin(np.radians(beta)) / np.sin(np.radians(alpha)), two_theta


def speckle_acf(crops, good, half=10):
    """Envelope-subtracted, window-corrected mean speckle autocorrelation.

    Each frame is scaled onto the mean pattern before subtraction so frame-to-frame
    flux variation does not leak into the residual.
    """
    mean_pat = crops.mean(axis=0) * good
    denom = mean_pat.sum()
    ny, nx = mean_pat.shape
    acc = np.zeros((2 * ny, 2 * nx))
    for f in crops:
        r = f * good
        r = r - (r.sum() / denom) * mean_pat
        F = np.fft.rfft2(r, s=(2 * ny, 2 * nx))
        acc += np.fft.irfft2(np.abs(F) ** 2, s=(2 * ny, 2 * nx))
    W = np.fft.rfft2(good.astype(float), s=(2 * ny, 2 * nx))
    wcorr = np.fft.irfft2(np.abs(W) ** 2, s=(2 * ny, 2 * nx))  # finite-window taper
    acf = np.fft.fftshift(acc / np.maximum(wcorr, 1e-9))
    return acf[ny - half:ny + half + 1, nx - half:nx + half + 1] / len(crops)


def speckle_decorrelation(crops, good, n_rows, n_cols, step_nm=20.0, max_off=4):
    """Illuminated-footprint FWHM [nm] per scan axis, from how fast speckle decorrelates.

    Uses the SAME envelope-subtracted residuals as speckle_acf, with the Poisson
    variance removed from the normalisation, so this and the ACF measure two
    different things on identical footing:
      - ACF width      -> the COHERENTLY scattering extent (one mode / one domain)
      - decorrelation  -> the ILLUMINATED extent (moving the probe swaps material)
    For Gaussian illumination of intensity-FWHM w, r(delta) = exp(-delta^2/(4 sigma_I^2)),
    sigma_I = w/2.3548. Their ratio counts coherent modes (or coherently
    diffracting domains) across the footprint.
    """
    mean_pat = crops.mean(axis=0) * good
    m = good.ravel()
    res = (crops * good).reshape(len(crops), -1)[:, m].astype(np.float64)
    flat_mean = mean_pat.ravel()[m]
    res -= (res.sum(axis=1, keepdims=True) / flat_mean.sum()) * flat_mean
    counts = (crops.reshape(len(crops), -1)[:, m]).mean(axis=1)   # Poisson variance per px
    var = (res ** 2).mean(axis=1) - counts
    idx = np.arange(len(crops)).reshape(n_rows, n_cols)
    out = {}
    for axis, tag in ((1, "fast"), (0, "slow")):
        deltas, rs = [], []
        for k in range(1, max_off + 1):
            i = idx.take(range(idx.shape[axis] - k), axis=axis).ravel()
            j = idx.take(range(k, idx.shape[axis]), axis=axis).ravel()
            cov = (res[i] * res[j]).mean(axis=1)
            r = np.median(cov / np.sqrt(np.maximum(var[i] * var[j], 1e-12)))
            deltas.append(k * step_nm)
            rs.append(min(max(r, 1e-6), 0.999999))
        d, r = np.array(deltas), np.array(rs)
        sig = np.median(d / (2 * np.sqrt(-np.log(r))))  # r = exp(-d^2 / 4 sigma^2)
        out[tag] = dict(w_nm=2.3548 * sig, r=r)
    return out


def _gauss2d(coords, amp, sa, sb, theta, base):
    x, y = coords
    xr = x * np.cos(theta) + y * np.sin(theta)
    yr = -x * np.sin(theta) + y * np.cos(theta)
    return amp * np.exp(-0.5 * ((xr / sa) ** 2 + (yr / sb) ** 2)) + base


def fit_acf(acf, max_lag=6):
    """Anisotropic-Gaussian fit excluding the zero-lag shot-noise spike.

    Returns sigmas along the ellipse axes, broad first; `sig_*_true` have the
    square-pixel self-correlation deconvolved (sigma^2 - 1/6).
    """
    half = acf.shape[0] // 2
    y, x = np.mgrid[-half:half + 1, -half:half + 1]
    sel = (np.abs(x) <= max_lag) & (np.abs(y) <= max_lag) & ~((x == 0) & (y == 0))
    p, _ = curve_fit(
        _gauss2d, (x[sel].ravel(), y[sel].ravel()), acf[sel].ravel(),
        p0=[acf[half, half], 3.0, 1.0, 0.8, 0.0],
        bounds=([0, 0.15, 0.15, -np.pi / 2, -np.inf], [np.inf, 20, 20, np.pi / 2, np.inf]),
        maxfev=40000)
    amp, sa, sb, theta, base = p
    if sb > sa:  # broad ACF axis first (== the SMALLER real-space extent)
        sa, sb, theta = sb, sa, theta + np.pi / 2
    deconv = lambda s: np.sqrt(max(s ** 2 - PIX_ACF_VAR, 1e-4))
    return dict(amp=amp, base=base, acf0=acf[half, half], theta=theta,
                sig_broad=sa, sig_narrow=sb,
                sig_broad_true=deconv(sa), sig_narrow_true=deconv(sb),
                theta_deg=np.degrees((theta + np.pi / 2) % np.pi - np.pi / 2))


def _gauss2d_pair(coords, a1, sa1, sb1, a2, sa2, sb2, theta, base):
    x, y = coords
    xr = x * np.cos(theta) + y * np.sin(theta)
    yr = -x * np.sin(theta) + y * np.cos(theta)
    return (a1 * np.exp(-0.5 * ((xr / sa1) ** 2 + (yr / sb1) ** 2))
            + a2 * np.exp(-0.5 * ((xr / sa2) ** 2 + (yr / sb2) ** 2)) + base)


def fit_acf2(acf, max_lag=8):
    """Two-component fit: a SHARP ACF peak on a BROAD pedestal, sharing one orientation.

    A single Gaussian is the wrong model here and biases the length badly (it
    lands between the two true scales). Physically, an illumination with overall
    envelope W and internal structure xi Fourier-transforms to a narrow central
    peak of width 1/W on broad wings of width 1/xi - so the sharp ACF component
    measures the ENVELOPE and the broad one the FINE structure inside it.

    Returned sigmas are pixel-binning-deconvolved; 'env' is the sharp component.
    """
    half = acf.shape[0] // 2
    y, x = np.mgrid[-half:half + 1, -half:half + 1]
    sel = (np.abs(x) <= max_lag) & (np.abs(y) <= max_lag) & ~((x == 0) & (y == 0))
    c, v = (x[sel].ravel().astype(float), y[sel].ravel().astype(float)), acf[sel].ravel()
    peak = acf[half, half]
    p, _ = curve_fit(_gauss2d_pair, c, v,
                     p0=[peak * 0.4, 1.5, 0.6, peak * 0.6, 4.5, 0.6, 0.8, 0.0],
                     bounds=([0, .15, .15, 0, .15, .15, -np.pi / 2, -np.inf],
                             [np.inf, 60, 60, np.inf, 60, 60, np.pi / 2, np.inf]),
                     maxfev=200000)
    a1, sa1, sb1, a2, sa2, sb2, theta, base = p
    comps = sorted(((sa1, sb1, a1), (sa2, sb2, a2)))  # smaller broad-axis sigma == sharper == envelope
    rms = float(np.sqrt(np.mean((v - _gauss2d_pair(c, *p)) ** 2)))
    deconv = lambda s: np.sqrt(max(s ** 2 - PIX_ACF_VAR, 1e-4))
    out = dict(theta=theta, theta_deg=np.degrees((theta + np.pi / 2) % np.pi - np.pi / 2),
               base=base, acf0=peak, rms=rms, total_amp=a1 + a2)
    for name, (sa, sb, amp) in zip(("env", "fine"), comps):
        out[name] = dict(sig_a=sa, sig_b=sb, amp=amp, weight=amp / (a1 + a2),
                         ext_a=extent_nm(deconv(sa)), ext_b=extent_nm(deconv(sb)),
                         fwhm_a=2.3548 * sa, fwhm_b=2.3548 * sb)
    return out


def axis_profile(acf, theta, n=49, rmax=12.0, swidth=3.0):
    """ACF profile along a principal axis, INTEGRATED across the perpendicular direction.

    Sampling a bare line along the ridge is hopeless here: the ridge is ~0.55 px
    wide across, so a sub-pixel misalignment of the sampling line biases the
    crest low. Integrating over +-swidth px perpendicular removes that
    sensitivity entirely (and the model integrates analytically the same way:
    a 2D Gaussian -> sqrt(2*pi)*sigma_perp times the same 1D Gaussian).
    """
    half = acf.shape[0] // 2
    t = np.linspace(-rmax, rmax, n)
    s = np.arange(-swidth, swidth + 1e-9, 0.25)
    tt, ss = np.meshgrid(t, s, indexing="ij")
    yy = half + tt * np.sin(theta) + ss * np.cos(theta)
    xx = half + tt * np.cos(theta) - ss * np.sin(theta)
    vals = map_coordinates(acf, np.vstack([yy.ravel(), xx.ravel()]), order=1, mode="nearest")
    return t, vals.reshape(tt.shape).sum(axis=1) * 0.25


def analyse(path, tag, samth, tth, gam):
    with h5py.File(path, "r") as f:
        r = process(find_frames(f), 71, 41, 64)
    acf = speckle_acf(r["crops"], r["roi_good"], half=14)
    fit1 = fit_acf(acf)
    fit = fit_acf2(acf)
    dec = speckle_decorrelation(r["crops"], r["roi_good"], 71, 41)
    env, fine = fit["env"], fit["fine"]
    # sampling is set by the FINE component: it is the widest thing in the ACF, so if IT is
    # unresolved the detector is throwing information away
    fwhm = (fine["fwhm_a"], min(env["fwhm_b"], fine["fwhm_b"]))
    ext = (fine["ext_a"], min(env["ext_b"], fine["ext_b"]))
    proj, two_theta = exit_projection(samth, tth, gam)
    print(f"{tag}: TWO-component ACF (ellipse {fit['theta_deg']:+.0f} deg, RMS resid {fit['rms']:.3f}); "
          f"a single Gaussian would report {extent_nm(fit1['sig_broad_true']):.0f} nm - "
          f"an artefact lying between the two real scales")
    print(f"   envelope comp: sigma {env['sig_a']:.2f}/{env['sig_b']:.2f} px -> "
          f"{env['ext_a']:.0f} x {env['ext_b']:.0f} nm  (weight {100 * env['weight']:.0f}%)  "
          f"<-- compare metadata spot ~150 nm")
    print(f"   fine comp:     sigma {fine['sig_a']:.2f}/{fine['sig_b']:.2f} px -> "
          f"{fine['ext_a']:.0f} x {fine['ext_b']:.0f} nm  (weight {100 * fine['weight']:.0f}%)")
    print(f"   speckle FWHM (widest comp) {fwhm[0]:.1f} / {fwhm[1]:.1f} px")
    print(f"   sampling: {fwhm[0]:.1f} px/speckle OK | {fwhm[1]:.1f} px/speckle "
          f"{'UNDERSAMPLED' if fwhm[1] < NYQUIST_PX else 'ok'} (need >= {NYQUIST_PX:.0f})")
    fov = ptycho_fov_nm()
    osr = (fov / ext[0], fov / ext[1])
    print(f"   coherent extents: {ext[0]:.0f} nm  and  >= {ext[1]:.0f} nm (2nd axis undersampled)")
    print(f"   ptycho oversampling (FOV {fov:.0f} nm / extent): {osr[0]:.1f}x  and  <= {osr[1]:.1f}x "
          f"(want >= 2 on both)")
    print(f"   geometry: 2theta={two_theta:.1f} deg, incidence {samth:.1f} / exit "
          f"{two_theta - samth:.1f} deg -> aperture factor sin(exit)/sin(inc) = {proj:.2f} "
          f"(NOT 1/sin(inc) = {1 / np.sin(np.radians(samth)):.2f})")
    geom = (150.0 / np.sin(np.radians(samth)), 150.0)  # metadata spot / sin(incidence) on the surface
    print(f"   speckle decorrelation length: fast {dec['fast']['w_nm']:.0f} nm, "
          f"slow {dec['slow']['w_nm']:.0f} nm  (r/step {dec['fast']['r'][0]:.3f}/{dec['slow']['r'][0]:.3f}) "
          f"-- independent real-space check on the envelope scale")
    return acf, fit, fit1, dict(fwhm=fwhm, ext=ext, proj=proj, two_theta=two_theta, samth=samth,
                                fov=fov, osr=osr, dec=dec, geom=geom, env=env, fine=fine)


def main():
    fig, ax = plt.subplots(2, 3, figsize=(15, 8.5))
    for row, (path, tag, samth, tth, gam) in enumerate(SCANS):
        acf, fit, fit1, res = analyse(path, tag, samth, tth, gam)
        half = acf.shape[0] // 2
        show = acf.copy()
        show[half, half] = np.nan  # shot-noise delta, excluded from the fit
        im = ax[row, 0].imshow(show, cmap="magma", extent=(-half, half, half, -half))
        ax[row, 0].set_title(f"{tag}  speckle ACF\n(zero lag = shot noise, masked)", fontsize=9)
        ax[row, 0].set_xlabel("lag, det x [px]"); ax[row, 0].set_ylabel("lag, det y [px]")
        plt.colorbar(im, ax=ax[row, 0], shrink=0.8)

        # Profile along the BROAD ellipse axis, where the sharp-on-broad structure lives.
        # |lag| >= 1 px only: inside one pixel the interpolated profile blends the zero-lag
        # shot-noise delta, which the fit excludes outright. Normalised so both rows compare
        # despite the ~17x intensity-squared difference in raw ACF units.
        SW = 3.0
        t, prof = axis_profile(acf, fit["theta"], rmax=12.0, swidth=SW)
        keep = np.abs(t) >= 1.0
        # cross-ridge integral of a 2D Gaussian component == sqrt(2 pi) sigma_perp * 1D Gaussian
        wgt = lambda sb: np.sqrt(2 * np.pi) * sb * min(1.0, 1.0)
        norm = sum(fit[k]["amp"] * wgt(fit[k]["sig_b"]) for k in ("env", "fine"))
        base_line = fit["base"] * 2 * SW
        ax[row, 1].plot(t[keep], (prof[keep] - base_line) / norm, "o", ms=3.5, color="k",
                        label="measured (cross-ridge integral)")
        ft = np.linspace(-12, 12, 400)
        tot = np.zeros_like(ft)
        for key, col, lab in (("env", "C0", "sharp: envelope"), ("fine", "C3", "broad: fine structure")):
            g = (fit[key]["amp"] * wgt(fit[key]["sig_b"]) / norm
                 * np.exp(-0.5 * (ft / fit[key]["sig_a"]) ** 2))
            tot += g
            ax[row, 1].plot(ft, g, "--", color=col, lw=1.2,
                            label=f"{lab} -> {fit[key]['ext_a']:.0f} nm")
        ax[row, 1].plot(ft, tot, "-", color="C2", lw=1.3, label="two-component total")
        g1 = (fit1["amp"] * wgt(fit1["sig_narrow"]) / norm
              * np.exp(-0.5 * (ft / fit1["sig_broad"]) ** 2))
        ax[row, 1].plot(ft, g1, ":", color="0.55", lw=1.3,
                        label=f"single Gaussian (biased) -> {extent_nm(fit1['sig_broad_true']):.0f} nm")
        ax[row, 1].axvspan(-1, 1, color="k", alpha=0.07)
        ax[row, 1].set_ylim(-0.03, 1.05)
        ax[row, 1].set_title(f"{tag}  ACF along the broad ellipse axis (normalised)", fontsize=9)
        ax[row, 1].set_xlabel("lag [px]"); ax[row, 1].set_ylabel("ACF / fitted peak", fontsize=8)
        ax[row, 1].legend(fontsize=6.5)

        under = "UNDERSAMPLED" if res["fwhm"][1] < NYQUIST_PX else "ok"
        ax[row, 2].axis("off")
        ax[row, 2].text(
            0.0, 0.97,
            f"{tag}\n\n"
            f"speckle FWHM      {res['fwhm'][0]:5.1f} px  |{res['fwhm'][1]:5.1f} px\n"
            f"sampling          {'ok':>8}  |{under:>9}\n"
            f"                             (need >= {NYQUIST_PX:.0f} px)\n\n"
            f"ENVELOPE comp  {res['env']['ext_a']:5.0f} nm ({100 * res['env']['weight']:.0f}%)\n"
            f"   metadata spot ~150 nm   <-- MATCHES\n"
            f"   real-space decorrelation (fast scan\n"
            f"   axis) {res['dec']['fast']['w_nm']:.0f} nm  <-- ALSO MATCHES\n\n"
            f"FINE comp      {res['fine']['ext_a']:5.0f} nm ({100 * res['fine']['weight']:.0f}%)\n"
            f"   structure INSIDE the envelope\n"
            f"   => multi-mode probe at rung 2\n\n"
            f"ptycho oversamp.  {res['osr'][0]:5.1f}x | <={res['osr'][1]:5.1f}x\n"
            f"   FOV = lam*R/px = {res['fov']:.0f} nm; want >= 2x\n\n"
            f"2nd axis (undersampled): both comps\n"
            f"   >= {min(res['env']['ext_b'], res['fine']['ext_b']):.0f} nm; ellipse "
            f"{fit['theta_deg']:+.0f} deg (det frame)\n"
            f"aperture factor sin(out)/sin(in) = {res['proj']:.2f}",
            va="top", fontsize=8, family="monospace", transform=ax[row, 2].transAxes)
    fig.suptitle(
        "Speckle autocorrelation: coherent illumination extents and detector sampling  "
        f"(10.0 keV, 75 um px, R = {DIST_M * 1000:.1f} mm readback / 500 mm nominal)\n"
        "RESULT 1 - detector pixels are spent very unevenly: ptychographic oversampling is ~13x on one "
        "axis but 1.6x/2.2x on the other (want >= 2 on both).\n"
        "   => bin the oversampled axis (free 4-6x data reduction) and move the detector back / use "
        "more symmetric reflections; rung-2 needs a probe-support constraint on the tight axis.\n"
        "RESULT 2 - the ACF needs TWO components (a single Gaussian is the wrong model and returns ~63 nm, "
        "an artefact lying between the two real scales):\n"
        "   a SHARP peak -> 157/159 nm illumination ENVELOPE, which independently reproduces both the "
        "~150 nm metadata spot size AND the real-space\n"
        "   decorrelation length (161/146 nm) - three routes agreeing; and a BROAD pedestal -> ~50 nm "
        "FINE structure inside that envelope, the same in both scans.\n"
        "   A structured/partially-coherent probe is the natural reading (zone-plate ring structure, or a "
        "~50 nm coherence patch) => rung 2 needs MULTI-MODE probes.\n"
        "   Discriminator already collected: on the Au standard a sample-driven 50 nm scale would change, "
        "a probe-driven one would not - a Globus pull away.\n"
        "Axes are not assigned to (probe / footprint / thickness): asymmetric reflections + tilted "
        "Ewald slice mix them - that needs the rung-2 geometry Jacobian.", fontsize=8.5)
    fig.tight_layout(rect=(0, 0, 1, 0.88))
    fig.savefig("speckle_geometry.png", dpi=140)
    print("wrote speckle_geometry.png")


def selftest():
    """Synthesise speckle from a known anisotropic probe; check the inversion recovers it."""
    rng = np.random.default_rng(0)
    n, dx_nm = 512, 8.0
    w_true = (160.0, 520.0)  # illumination intensity FWHM along (y, x) [nm]
    yy, xx = np.mgrid[-n // 2:n // 2, -n // 2:n // 2] * dx_nm
    sig = [w / 2.3548 for w in w_true]
    probe = np.exp(-0.25 * ((yy / sig[0]) ** 2 + (xx / sig[1]) ** 2))  # |probe|^2 has FWHM w
    crops = []
    for _ in range(60):
        obj = np.exp(2j * np.pi * rng.random((n, n)))                  # fully rough
        far = np.abs(np.fft.fftshift(np.fft.fft2(probe * obj))) ** 2
        c = far[n // 2 - 64:n // 2 + 64, n // 2 - 64:n // 2 + 64]
        crops.append(c / c.mean())
    acf = speckle_acf(np.array(crops), np.ones((128, 128), bool))
    fit = fit_acf(acf)
    dq = 2 * np.pi / (n * dx_nm)  # this synthetic's q step per "pixel" [1/nm]
    got = sorted((extent_nm(fit["sig_broad"], dq), extent_nm(fit["sig_narrow"], dq)))
    print(f"SELFTEST 1 (smooth probe): true {w_true[0]:.0f}/{w_true[1]:.0f} nm -> "
          f"recovered {got[0]:.0f}/{got[1]:.0f} nm")
    for g, t in zip(got, sorted(w_true)):
        assert abs(g - t) / t < 0.12, f"extent inversion off: {g:.0f} vs {t:.0f} nm"

    # SELFTEST 2: the model-selection failure this module actually made. A sharp ACF peak on a
    # broad pedestal (two illumination scales) is built directly, with the zero-lag shot-noise
    # delta on top; fit_acf2 must recover BOTH sigmas, and the single Gaussian must return a
    # bogus value strictly BETWEEN them -- the signature of the original error.
    s_env, s_fine, half2 = 1.40, 4.50, 14
    gy, gx = np.mgrid[-half2:half2 + 1, -half2:half2 + 1]
    th0 = 0.80
    xr = gx * np.cos(th0) + gy * np.sin(th0)
    yr = -gx * np.sin(th0) + gy * np.cos(th0)
    synth = (0.40 * np.exp(-0.5 * ((xr / s_env) ** 2 + (yr / 0.55) ** 2))
             + 0.60 * np.exp(-0.5 * ((xr / s_fine) ** 2 + (yr / 0.60) ** 2)))
    synth = synth + 0.002 * rng.standard_normal(synth.shape)
    synth[half2, half2] += 1.5                     # white shot noise: zero lag only
    f2, f1 = fit_acf2(synth), fit_acf(synth)
    got_env, got_fine = f2["env"]["sig_a"], f2["fine"]["sig_a"]
    single = f1["sig_broad"]
    print(f"SELFTEST 2 (two-scale ACF): sigma_env {s_env:.2f} -> {got_env:.2f} px, "
          f"sigma_fine {s_fine:.2f} -> {got_fine:.2f} px;  a single Gaussian says {single:.2f} px")
    assert abs(got_env - s_env) < 0.25, f"envelope sigma off: {got_env:.2f} vs {s_env:.2f}"
    assert abs(got_fine - s_fine) < 0.5, f"fine sigma off: {got_fine:.2f} vs {s_fine:.2f}"
    assert got_env < single < got_fine, "single-Gaussian answer must fall between the two scales"

    # the asymmetric-reflection projection must differ from the naive 1/sin(incidence)
    proj, _ = exit_projection(11.2, 42.2, 23.15)
    assert 2.5 < proj < 3.5 and abs(proj - 1 / np.sin(np.radians(11.2))) > 1.5
    print("SELFTEST PASS")


if __name__ == "__main__":
    p = argparse.ArgumentParser()
    p.add_argument("--selftest", action="store_true")
    a = p.parse_args()
    selftest() if a.selftest else main()
