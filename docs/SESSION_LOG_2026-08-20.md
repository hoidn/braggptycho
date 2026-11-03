# Session log — NbN Bragg-nanodiffraction problem-space analysis

**Date:** 2026-08-20  **Repo:** `/home/ollie/Documents/braggptycho`
**Collaboration:** SLAC (Oliver) + ANL, APS beamline 26-ID; epitaxial δ-NbN(111) superconducting thin films.

## 0. The task

> Read ctx.md. Do a scientific analysis of the problem space. Understand the scientific intent of the
> people involved and the essence of the experimental, data-collection, imaging, and analysis
> challenges. Then brainstorm multiple possible directions and suggest a couple of options for an
> executable plan that would generate progress to show at the upcoming meeting. Retrieve and review
> datasets and relevant papers, if possible.

Follow-ups during the session: (a) get a subagent onto the *actual* data; (b) add a dedicated
S. O. Hruszkewycz publication survey for literature context; (c) surface/send the "money figure";
(d) make `twin_comparison.png` a permanent repo artifact; (e) this log.

## 1. Deliverables produced

| File | What it is |
|---|---|
| `BRIEFING.md` | The full written deliverable: science + each collaborator's goals (§1), measurement/imaging challenges (§2), real-data findings (§3), the staged direction ladder (§4), the two meeting plans (§5), asks for David (§6), references (§7). |
| `quicklook.py` | Tested reconstruction-free first-look pipeline. Masks the bad-pixel sentinel, infers the raster grid, computes mean pattern / ROI / peak-COM / integrated-intensity / adjacent-frame speckle-correlation maps. Has a synthetic-twin `--selftest` that passes. |
| `compare_twins.py` | Renders the meeting "money figure" `twin_comparison.png` from the two `*_maps.npz`. |
| `twin_comparison.png` | The money figure (2×3 panels: intensity, COM_y, speckle-corr for (202) vs (131)). Now permanent in-repo. |
| memory: `braggptycho-project.md` + `MEMORY.md` | Durable project memory (goals, participants, verified data facts, direction ladder). |
| scratchpad: `dossier.md`, `empirical_findings.md`, `data/scan_*_maps.npz`, `data/scan_*_quicklook.png` | Fact bases + reduced maps + per-scan quicklook figures. Source `.h5` (≈700 MB each) live only in scratchpad. |

## 2. What was done, in order

1. **Read `ctx.md`** (ANL meeting 1 notes, SLAC meeting 2 auto-transcript, David Garcia's data email).
   Established the collaboration structure, science goals, scan numbers (308/309 = (202), 312/313 =
   (131)), the h5 layout, and the Globus/OneDrive access route.
2. **Ran the `deepthink` 14-step workflow** to structure the problem-space analysis; carried it to
   CERTAIN confidence.
3. **Launched 8 subagents:** 3 literature (general BPP/CDI, reconstruction-free methods, dose/info
   theory), 1 data-hunt, and 4 analysis personae (physicist, algorithmist, strategist, skeptic). Later
   added a **dedicated Hruszkewycz survey** agent (completed post-compaction; see §6).
4. **Retrieved 2 of 4 real `.h5` files** (scans 309 and 313, the 10 nm film) via the data-hunt agent,
   which used the user's *live* Firefox SharePoint session. Read-only; the cookie artifacts were
   deleted immediately afterward.
5. **Built and tested `quicklook.py`.** The self-test synthesises a twin-domain scan (bulk lattice
   tilt + fine-scale disorder for genuine speckle) and asserts grid recovery, correlation decay, and a
   domain-COM split. Passes.
6. **Ran on the real data**, iterated through several data-reduction bugs (§4), and produced the
   verified reconstruction-free maps + the money figure.
7. **Wrote `BRIEFING.md`** with the two executable plans and the asks for David.
8. **Hruszkewycz survey returned** — corrected two mis-attributions, confirmed "Steve," and surfaced
   two high-value papers now folded into `BRIEFING.md §7`.

## 3. Key decisions (and the reasoning)

- **Stay in the 2D single-angle regime for thin films.** For a film of thickness *t*, the Bragg peak
  *is* the q_z rod (width ≈ 2π/t); depth-resolved single-angle 3D BPP is fragile and unnecessary here.
  Hruszkewycz's own Nat Mater 2017 is cited as the *fragility* argument, not as something to emulate.
- **Reconstruction-free maps first (Plan A leads), learned reconstruction as a parallel seed (Plan B).**
  Plan A delivers real physics on the shared data *now* with no positions file; Plan B (PtychoPINN-Bragg
  + Fisher-information acquisition) is the user's genuine research niche but is blocked on metadata, so
  it seeds in parallel rather than gating the meeting.
- **Mask the bad-pixel sentinel before any arithmetic.** Non-negotiable: unmasked sums are ~8 orders of
  magnitude wrong (§4).
- **Resolve grid orientation by map-coherence, not by assumption.** The auto-resolver picks the
  orientation that minimises total-map roughness. This *reversed* an early wrong guess (§4).
- **Present COM as a *mixed* tilt+strain map, explicitly not a clean strain map.** A single fixed-θ scan
  gives 2 numbers/pixel and cannot separate tilt (2 DOF) from strain (1 DOF). Honest labelling avoids
  overclaiming; the clean split needs a matched-raster multi-angle re-scan (Holt PRL 2014 decomposition).
- **Do NOT silently resolve the energy discrepancy — flag it to David.** The file's `DCM Energy` field
  (12.5 keV) is inconsistent with the diffraction geometry, which indexes cleanly only at 10.0 keV. This
  is surfaced as an open question, not patched over (specs > implementation).
- **Treat (202) and (131) as different *regions* (radiation damage forced region-hopping)** → only a
  statistical comparison between them is valid, never per-pixel; cross-peak u(r) fusion is impossible on
  *this* dataset and is deferred to a future matched-raster beamtime.
- **Delete the SharePoint cookie artifacts after the read-only data pull.** Minimise credential
  footprint; the user's email/identity is for identification only.

## 4. Data-reduction bugs found and fixed (real cost centres)

- **Bad-pixel sentinel = 4294967295 (2³²−1)**, 3.58 % of pixels including whole detector module-gap
  rows. Unmasked → results ~8 orders of magnitude wrong. Fix: mask first, everywhere.
- **Raster grid is 71×41, not 41×71.** An early forced `--grid 41,71` was *transposed* and manufactured
  fake horizontal stripes. Caught by a serpentine/roughness test (roughness 0.449 for 71×41 vs 1.226 for
  41×71). The auto-resolver now picks 71×41 once masking is applied.
- **No serpentine artifact** (even/odd row COM differ < 0.05 px; serpentine "correction" makes maps
  worse — so it's off).
- **No per-frame positions in the shared files** (ZPX/ZPY are scan-centre scalars; Positioner absent;
  NDAttributes empty). This *blocks* position-indexed reconstruction → a concrete ask for David.
- **Beam-center argmax locked onto a flagged pixel** (row 49) until masking was applied; the real (202)
  peak is at (813, 591).
- **Filename bug** (`1003${s}` → `scan_309_1003309.h5`) fixed to the correct `scan_309_100309.h5`.
- **Self-test tuning:** enlarged synthetic domains (~300 nm) + domain-correlated bulk tilt so the COM
  split is visible; added ~30 nm Gaussian disorder for genuine speckle decay; switched the decay
  assertion to an aggregate (element-wise monotonicity legitimately fails at domain-boundary tilt flips
  — a *real* physical effect, not a bug).

## 5. Empirical findings on the real data (scans 309/313, 10 nm film)

- **(202) scan 309:** bright, integrated ROI ~260k ph/frame. Integrated intensity is **flat**
  (max/median ≈ 1.1) but spatially **coherent** and heavily speckled (adjacent-frame corr ≈ 0.90,
  anisotropic: fast 0.916 / slow 0.905) — matches David's "flat intensity, lots of speckle." COM maps
  show coherent lattice tilt/strain domains, hundreds-of-nm, with sharp boundaries.
- **(131) scan 313:** ~4× dimmer (integrated ~58–73k). Integrated intensity is **strongly structured** —
  diagonal/wavy bright–dark bands (the twin/stripe morphology). Speckle correlation is **lower** (fast
  0.689 / slow 0.677) — a direct, visible signature of **Poisson photon-starvation** decorrelating the
  weak reflection.
- **Geometry (physicist-verified):** footprint elongation 5.15× (202) / 3.14× (131); sample-frame grid =
  zpx/sin θ → ~5.1× x-pitch anisotropy at θ = 11.2°; q_z rod width 1.05/0.63/0.16 nm⁻¹ for 6/10/40 nm →
  specular rocking FWHM 1.22/0.73/0.18° (0.3° step is marginal at 10 nm, ~2.4 samples/FWHM); real-space
  resolution ~15–20 nm from the detected q-extent; ZP2 badly out of focus at large FOV.
- **Damage mechanism:** NbN absorbs only ~1.5×10⁻⁴ of incident photons → damage is an event-cascade
  (N-vacancy) process, not bulk heating; empirical single-region ceiling ~1 h 50 min (Night5Ptycho uses
  ~39 min/region, ~2× margin) → dose budgeting, not thermal.

## 6. Literature — the S. O. Hruszkewycz survey (the last piece)

Confirmed **"Steve" ≈ Stephan O. Hruszkewycz** (ANL MSD group leader, inventor of nanofocused BPP,
works at 26-ID since 2011) — inference, high confidence; competing ANL "Steves" don't fit.

**Two attribution corrections** (aimed at the raw digests; `BRIEFING.md` already had Burdet right):
- The "BPP on a Nb film, PRB 2017" precedent is **Burdet et al. (Robinson group), PRB 96, 014109 (2017)**
  — *not* Hruszkewycz. Still the nearest experimental precedent; cite it, attribute it correctly.
- The Annu. Rev. Mater. Res. 2013 nanoprobe tutorial is **Holt et al.**, not him.

**Two high-value papers I had missed — now added to §7:**
- **Hruszkewycz PRL 110, 177601 (2013):** BPP phase → *picometre displacement / polarization map* of
  180° nanodomains in an epitaxial film. This is the **direct template** for our twin-domain
  displacement deliverable — the closest published analog to the project's central goal.
- **Calvo-Almazán Sci. Rep. 9, 6386 (2019):** jointly reconstructs the crystal *and each frame's true
  angle* — i.e. a **software answer to our missing/uncertain scan-position problem**, the exact blocker
  we hit on the shared data.

Also folded in: **Zhou PRB 2024** (Takagi–Taupin for focused-nanobeam Bragg — flags when kinematic BPP
phase is biased by the dynamically-diffracting substrate under an epitaxial film), the corrected
**Kandel Opt. Express 32, 19594 (2024)** citation for multi-peak BPP, and **Porter JAC 2026** (newest
robust multi-peak coupling).

**Positioning insight (now in §7):** his corpus has *every* ingredient — nanofocused BPP on
thin films, phase→displacement, defect/stacking-fault phasing, angular-registration recovery, multi-peak
vector imaging — but **no NbN or superconducting-film paper**. This project would be his program's first
entry into superconducting quantum-device materials: a clean novelty argument for the proposal.

## 7. The direction ladder (one pipeline, staged — not competing tracks)

0. **Reconstruction-free maps** — *DONE* (intensity / peak-COM / speckle-corr, both twins).
1. **COM → gradient-integrate → smoothed G·u** (defuses phase-ramp pile-up).
2. **Single-peak 2D BPP** — no off-the-shelf engine; build a geometry layer in front of pty-chi /
   PtychoPINN.
3. **PtychoPINN-Bragg** (the open niche: overlap + Poisson loss vs. the per-pattern DONUT/NanobeamNN
   approach) + **Fisher-information acquisition** (optimal rocking offset Δ* = σ√2 off-peak).
4. **Multi-peak fusion** — deferred to the next beamtime (needs matched-raster multi-angle scans).

## 8. Open questions / asks for David

- **Scan positions / master file** (per-frame ZPX/ZPY) — blocks any reconstruction.
- **Pixel size + beam center** — absent from the shared files.
- **Confirm the operative energy** (geometry says 10.0 keV; the file's `DCM Energy` says 12.5 keV).
- For a future beamtime: **matched-raster multi-angle rocking** on the *same* region, to (a) separate
  tilt from strain cleanly and (b) enable cross-peak u(r) fusion.

## 9. Standing caveats (state these at the meeting)

- Grid 71×41 is *inferred* (map-coherence); step/pixel units are approximate pending David's metadata.
- (202) vs (131) are **different regions** (damage) → statistical comparison only, never per-pixel.
- COM conflates mosaic tilt + d-spacing strain → labelled "mixed tilt-dominant," not a clean strain map.
- `quicklook.py` has benign pyright warnings from the h5py stubs (safe to ignore).
- `compare_twins.py` reads the `*_maps.npz`, which currently live only in scratchpad; copy them into
  `data/` if the figure must be regenerable from the repo alone (the 700 MB `.h5` stay in scratchpad).

## 10. Addendum — adversarial review + repair (same day, second session)

A decision-critic pass verified every checkable claim against the raw `.h5` and David's experiment
notebooks (scan sheet CSV + log doc, in the first session's scratchpad). Outcomes:

**Confirmed / strengthened**
- Energy settled at **10.00 keV**: scan-sheet angles (309: samth 11.2°, tth 42.2°, γ 23.15°;
  313: samth 18.6°, tth 52.4°) index (202) and (131) to a=4.39 Å NbN at exactly 10.00 keV
  (cos2θ = cos tth · cos γ). `DCM Energy` 12.5 keV is stale metadata.
- Overnight cadence verified: scans 308–313 start ~39 min apart, zpy stepped 10 µm/scan (fresh region).
- **Poisson-starvation claim upgraded from plausible to proven**: shot-noise-corrected Pearson
  (subtract Poisson variance = per-frame mean from the normalization) sends both reflections to the
  same underlying correlation — (202) 0.998 fast/0.987 slow, (131) 0.999/0.983 vs raw 0.90/0.89 and
  0.66/0.65 — with raw values predicted to ~1% by the signal-variance fraction. Now computed by
  `pearson_maps()` (`corrc_*`, `decayc_*` in the npz) and shown in the twin-comparison figure.

**Corrected (previous claims wrong or overstated)**
- **The "~5.15× sample-frame x-pitch anisotropy" was most likely backwards.** The beamline macros in
  the scan sheet pre-compress zpx by sinθ (`-30*sin(samth)..`), making the *sample-frame* grid
  isotropic; the data agree (per-step decorrelation ~5× *smaller* along fast/x — an elongated
  *footprint*, not an elongated pitch; the two cancel in fractional-overlap terms, so equal-pitch
  would predict equal decorrelation). What is anisotropic is the probe **footprint/resolution**
  (~150 nm y × ~0.5–0.8 µm x). Pitch stays "nominal 20 nm" until David confirms `Night5Ptycho()`.
- "(131) strongly structured vs (202) flat" softened to numbers: std/med 6.5% vs 3.0% (~2×).
- The first-session npz were made with undocumented `--roi-half 80`; maps regenerated into repo
  `data/` with default flags, and all quoted numbers updated to the reproducible values
  (229k/55k counts/frame; 14/3.4 ph/px).
- Index-vs-position correlation concern (raised by Oliver): checked — `pearson_maps` pairs frames only
  within grid rows/columns after the reshape (no row wraparound); validity rests on the raster
  assumption, which is triple-checked (roughness, serpentine test, COM coherence) and fully resolved
  only by the positions file (top ask). Do NOT 1/sinθ-rescale the decay curves — that would bake in
  the probably-wrong pitch reading.

**Built (previously promised, didn't exist)**
- `fisher_rocking.png` (Bet B2 seed: σ√2 optimum) and `photon_budget.png` (Bet B1 seed: radial
  counts/pixel, 10 ph/px by r≈50 px, 1 ph/px by r≈105 px).
- `twin_comparison.png` overhauled to 2×4 (intensity, COM, raw corr, corrected corr) with visible
  correlation scales (0.5–1.0) and a corrected title (no pitch-anisotropy assertion; logged angles
  + 10.0 keV added).

## 11. Addendum 2 — `screenshots/` pass (David's meeting slides + notebook captures)

14 screenshots (2026-08-07, the SLAC meeting) were reviewed against our claims. Net effect:

- **Pitch convention upgraded to confirmed:** the sheet's fly2d rows all use
  `zpx, ±X*sin(radians(samth.VAL))`, including **scan 289 on the 10 nm sample at the (131) angle**
  (zpx ±1.914 = ±6·sin 18.6°). Sample-frame pitch is isotropic; only Night5Ptycho's one-line confirm
  remains. Figure/BRIEFING/memory updated.
- **Reconstruction reclassified soft-blocked:** the beamline "Diffraction User Data Explorer" table
  records M1/M2 MIN/MAX/NP for every scan (incl. macro-internal). Rows 308–313 + the convention give
  exact nominal position grids; Calvo-Almazán-style refinement absorbs the residual. Ask list updated
  (screenshot of those rows beats a master-file transfer).
- **Attribution corrections (avoid reinventing/stating the obvious at the meeting):** the (131)
  stripe-domain intensity morphology is David's prior observation (scan 289 "is this real?? yes",
  fly maps 305–307) — our map is a cross-validation at 20 nm steps, not a discovery. Virtual-aperture
  darkfield maps are his live GUI workflow (his satellite maps) — dropped as a "new" Plan-A item.
  The shot-noise correlation correction is standard XPCS practice — cite it as such; the new part is
  using it to separate structure from starvation across twin reflections.
- **New ensemble context folded into §1:** substrate is sapphire (twins nucleate on its O sublattice);
  6 nm film shows opposite in-plane twists (~176°/184°) locking to 180° by 40 nm; GIXRD-3D RSM
  diffuse-scattering analysis exists — our real-space maps extend his narrative, not replace it.
- **Confirmed details:** 0.5 s/frame exposure; detector `s26_eiger2_cnm` (Eiger2 → 75 µm px);
  detector distance still ambiguous (522.5 vs 395/400/500 mm notes) — kept as ask.
- New idea logged (optional): per-frame speckle-orientation (structure-tensor)
  map of the in-peak streak texture = local fault-plane orientation — the PI's literal "rotate so the
  lines are parallel" idea; no beamline tool does it.

## 12. Addendum 3 — direction reevaluation + rung-1 execution (`gu_proxy.py`)

With positions soft-unblocked, rung 1 (COM → integrate → G·u) became executable and was executed.
**The data revised the rung's own premise.** Journey (kept honest for the record):

- Naive 8-pairing integration: 80–85% of COM-field power non-integrable.
- Free 2×2 detector→sample map fit (generalized eigenproblem): degenerate — collapses onto the stripe
  direction (rank-1), rejected.
- Rotation-only frame fit + 2-step low-pass: residual barely moves (smoothing does NOT remove the
  curl) → the non-integrability is **domain-scale structure, not speckle noise**.
- Final form: Helmholtz decomposition at the most-gradient-compatible frame. **≤44% ((202)) / ≤34%
  ((131)) of domain-scale COM power can be a scalar strain-gradient surface; the majority is
  rotational — real-space twist-dominated mosaic**, matching David's ensemble twist-mosaicity story.
  Deliverables: twist-domain map Ψ (new modality) + gradient-part G·u proxy (±25 pm, provisional
  calibration: Eiger2 75 µm @ R=522.5 mm, 10 keV, 20 nm steps). Consistent bonus: corr(twist,
  intensity) = −0.24/−0.25 in both scans (independent regions/reflections) — suggestive only.
- Physics notes for interpretation: per-frame COM in the fully-speckled regime = envelope tilt +
  speckle-realization noise; each reflection sees only its own twin, so COM in the other twin's
  domains is dark-region noise; piecewise tilt fields carry curl on domain walls — all three push
  against scalar-phase (DPC) shortcuts and toward phasing. BRIEFING §4 rung 1 rewritten accordingly.
- Not built (deliberately): actual rung-2 reconstruction attempt (high failure risk at this stage) and the
  structure-tensor orientation map (twist map fills the new-modality slot; orientation map stays
  queued). Scans 308/312 (40 nm-step twins, in David's OneDrive links) would add a second step-size
  point to the correlation-vs-overlap information story — optional future pull.

## 13. Addendum 4 — speckle geometry (`speckle_geometry.py`), and two claims it corrected

Triggered by the notebook finding that the detector distance is ambiguous (R readback 522.5 = nominal
500 mm + ~22.8 mm), which makes every q/pm calibration provisional. Rather than wait on David, the
speckle autocorrelation measures the illuminated geometry from the data itself. Inversion
(`w = 1.665/(sigma_ACF · dq_per_px)`, Siegert) is validated by `--selftest`: synthetic speckle from a
known 160/520 nm anisotropic probe is recovered as 159/522 nm.

**Results**
- **Detector sampling (interpretation-independent):** speckle FWHM 8.5 px on one axis,
  1.4–1.7 px on the other — below the 2 px Nyquist floor. Ptychographic oversampling (FOV = λR/px =
  864 nm) is 13.5× vs 1.6× ((202)) / 2.2× ((131)). Actions: probe-support constraint at rung 2, free
  4–6× binning of the oversampled axis, detector further back next beamtime.
- **Two illumination length scales; probe size confirmed.** *(This bullet supersedes a wrong first
  version — see §14.)* The ACF needs two components: a sharp peak → **157/159 nm envelope**, and a
  broad pedestal → **~50 nm fine structure inside it** (weights ≈42/58, 46/54; two-component fit cuts
  the residual 29%/50%). The envelope is confirmed three ways — ACF sharp component, David's ~150 nm
  metadata, and real-space decorrelation vs probe translation (161/146 nm). ⇒ structured or partially
  coherent probe ⇒ **multi-mode probes at rung 2**. Discriminator already collected: a sample-driven
  50 nm scale must change on the Au standard and across the 6/40 nm films; a probe-driven one will not.

**Two earlier claims corrected by this pass**
- The far-field aperture does **not** scale as 1/sin(incidence). These are asymmetric reflections
  (2θ = 47.1°, 11.2° in / 35.9° out; and 55.9°, 18.6°/37.3°), so the footprint is projected onto the
  plane normal to the *exit* beam: factor sin(exit)/sin(incidence) = **3.02 / 1.90**, not 5.15 / 3.14.
  The 1/sin(incidence) elongation remains correct for the *surface footprint* and hence map resolution
  — the two were being conflated. BRIEFING §2 now separates them explicitly.
- The isotropic-pitch evidence was quoted as "per-step decorrelation ~5× smaller along fast/x". The
  measured ratio is **1.7–2.0×**, not 5. The argument still holds (equal-pitch-with-elongated-footprint
  predicts anisotropy; the "5× coarser x-pitch" reading predicts *none*, and the ordering across the
  two reflections is right), but the magnitude claim was wrong and is fixed.

## 14. Addendum 5 — a wrong model caught by inspection, and the guard added

Oliver looked at the ACF profile panel and asked why there was "a sharp peak on top of a broader one."
I had dismissed those points as an interpolation artefact. He was right and I was wrong.

**The error.** The speckle ACF was fitted with a *single* anisotropic Gaussian. The real ACF has two
components, so the single fit landed *between* the two true scales and returned ~63 nm — from which I
concluded "the coherently diffracting region is 2–5× smaller than the illuminated footprint," offering
partial coherence vs small mosaic domains as explanations. **The gap was an artefact of the model, not
a property of the sample.**

**The correct result.** Two-component fit (RMS residual −29% for (202), −50% for (131)), with both
scans independently agreeing: sharp peak → 157/159 nm envelope; broad pedestal → 50/49 nm fine
structure. Physically standard: an illumination with envelope *W* and internal structure *ξ* has an ACF
that is a narrow peak of width 1/*W* on broad wings of width 1/*ξ*. The envelope now agrees with
David's ~150 nm metadata *and* with the independent real-space decorrelation length (161/146 nm) —
three routes to one number, where before there was an unexplained factor of 2–5.

**Why the two nuisance explanations are excluded:** the sharp component is anisotropic (σ 1.5 × 0.55 px),
so it is not isotropic detector charge sharing; and it has width well beyond the single zero-lag pixel,
so it is not the white shot-noise delta (which the fit excludes anyway).

**Guard added.** `--selftest` gained a second case that builds a two-scale ACF with a zero-lag noise
spike and asserts that (a) `fit_acf2` recovers both σ (1.40→1.41, 4.50→4.52 px) and (b) the single
Gaussian returns a value *strictly between* them (3.69 px) — i.e. the failure mode is now reproduced
and asserted against. A first attempt at this selftest was itself wrong (the "fine structure" I
synthesised had a ~400 nm correlation length, so there were no two scales to separate); rebuilt to
test the fitting machinery directly, which is what had failed.

**Figure fixes from the same pass:** the profile is now the *cross-ridge integral* rather than a bare
interpolated line (the ridge is 0.55 px wide across, so sub-pixel misalignment of a sampling line
biases the crest low — that, not the noise delta, was distorting the plotted core), and the panel no
longer prints the retracted "coherent region smaller" line.

**Discarded during the pass (recorded so it isn't retried):** a first version read the ACF axes
directly as "probe = 63 nm" and tested the aspect ratio against 1/sinθ. Rejected — the narrow axis is
sub-pixel (so its aspect is extrapolation), the projection factor is the asymmetric one, and the
tilted Ewald slice mixes probe/footprint/thickness onto the detector axes. Axis identification is
deferred to the rung-2 Jacobian; the figure states this. An intermediate "~3 coherent modes per
footprint" line was also cut: for partially coherent illumination the ACF width and the translation
decorrelation both track the coherence patch, so their ratio is not a mode count.
