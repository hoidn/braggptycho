#!/usr/bin/env python3
"""Plan-B seed: where on the rocking curve does one photon buy the most strain info?

For a Gaussian rocking curve I(th) = I0 exp(-th^2/2s^2) and Poisson counting,
the Fisher information about a peak shift mu (== strain along G) from measuring
at fixed offset th is F(th) = (dI/dmu)^2 / I = I0 (th/s^2)^2 exp(-th^2/2s^2):
zero AT the peak, maximal at th = s*sqrt(2). Answers meeting-1's literal
question "where on the rocking curve do you want to be?" (step size / overlap
held fixed; per incident photon, kinematic, single peak).
"""
import numpy as np
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt

th = np.linspace(-3.5, 3.5, 1401)  # offset in units of rocking-curve sigma
I = np.exp(-th ** 2 / 2)
F = th ** 2 * I  # Fisher info per unit incident flux (arb. units), sigma=1
opt = np.sqrt(2)
assert abs(abs(th[np.argmax(F)]) - opt) < 5e-3, "Fisher max must sit at sigma*sqrt(2)"

fig, ax = plt.subplots(figsize=(8, 5))
ax.plot(th, I, "C0", lw=2, label="rocking curve  I(th)/I0")
ax.plot(th, F / F.max(), "C3", lw=2, label="strain Fisher info / photon (norm.)")
ax.axvline(0, color="gray", lw=0.8)
for s in (-opt, opt):
    ax.axvline(s, color="C3", ls="--", lw=1)
ax.annotate("zero information\nexactly at the peak", xy=(0, 0.02), xytext=(0.25, 0.35),
            arrowprops=dict(arrowstyle="->"), fontsize=10)
ax.annotate(r"optimum $\Delta^*=\sigma\sqrt{2}\approx1.41\sigma$", xy=(opt, 1.0),
            xytext=(1.7, 0.8), arrowprops=dict(arrowstyle="->"), fontsize=10, color="C3")
ax.set_xlabel(r"rocking offset  $\Delta\theta$  [$\sigma$ of rocking curve]")
ax.set_ylabel("normalized")
ax.set_title("Dose-limited acquisition: measure OFF the rocking-curve peak\n"
             "(Gaussian curve, Poisson counting; overlap/step fixed; CRB on strain = 1/sqrt(N F))",
             fontsize=11)
ax.legend(loc="upper left", fontsize=9)
fig.tight_layout()
fig.savefig("figures/fisher_rocking.png", dpi=140)
print("wrote figures/fisher_rocking.png; argmax at", round(float(abs(th[np.argmax(F)])), 3), "sigma")
