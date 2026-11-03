# NbN Bragg-nanodiffraction collaboration — analysis & plan

**Prepared 2026-08-20; reviewed and revised the same day.**
Working files in this repo: `quicklook.py` (tested pipeline, `--selftest` passes),
`compare_twins.py` → `twin_comparison.png`, `fisher_rocking.py` /
`photon_budget.py` → the two Plan-B seed figures, reduced maps in `data/*_maps.npz`
(regenerated with default flags — quoted numbers reproduce). Two real files pulled from
David's OneDrive: `scan_309` = NbN (202), `scan_313` = NbN (131), 10 nm film; raw `.h5`
(~700 MB each) are not in the repo. David's experiment notebooks (scan sheet + log)
were also retrieved and are now the source for the verified geometry below.

---

## 1. The science, and what everyone wants out of it

**System.** Epitaxial δ-NbN(111) superconducting films **on sapphire**, thickness series
6 / ~10-12 / 40 nm, straddling a layer-by-layer → columnar growth-mode transition near 10 nm.
NbN(111) on the sapphire **oxygen surface sublattice** carries a fundamental **180° twin degeneracy**
(fcc ABC vs ACB stacking); the two variants have equal bonding energy, so twins nucleate densely and
stochastically. David's GIXRD-3D RSM slides (in `screenshots/`) frame the ensemble picture: in the
6 nm film the two variants are **twisted in opposite in-plane directions** (~176°/184°, "twist
mosaicity"), locking to a clean 180° by 40 nm while a monoclinic distortion grows; 25-nm satellite
ordering in the 6 nm film, few-hundred-nm stripe domains by 12 nm, and ordered-defect
("kaleidoscopic") diffuse scattering — stacking-fault rods, partial dislocations — in the 40 nm film.
He presents this as strain accumulating in twinned domains, then relieving into defects past the
critical thickness. The nanoprobe work is the real-space, spatially resolved version of that story and
builds on his ensemble picture.

**Why it matters.** These are exactly the defects that set Tc, order-parameter uniformity, and device
yield in the two live NbN technologies — SNSPDs (5-10 nm films, i.e. this thickness window) and
all-nitride superconducting qubits. The field currently treats "disorder" as a scalar (sheet
resistance); nobody has spatially mapped the twin/strain/fault populations on the as-grown film at the
relevant length scale. That is the scientific opening.

**What each collaborator wants:**
- **David (grower):** the data gets looked at and shows recognizable structure vs thickness.
- **Senior PI:** his "4 numbers per image" correlation-map idea, tried on real data — does a
  reconstruction-free map resemble real space?
- **ANL Bragg-ptychography experts ("Steve" is very plausibly Stephan Hruszkewycz):** the
  reconstruction path is unblocked; eventually a displacement field u(r) reconciling the peaks.
- **Oliver (ML/algorithms):** PtychoPINN and the Fisher-information framing as a distinct
  contribution.

## 2. The measurement, and why it's a hard imaging problem

Hard X-ray nanoprobe at APS 26-ID, zone-plate raster (sample fixed). The essential structure:

- **Thin film ⇒ 2D is the right regime.** Each Bragg peak is a q_z rod of width 2π/t (t=6/10/40 nm →
  1.05/0.63/0.16 nm⁻¹). Single-angle **2D Bragg projection ptychography (BPP)** is well-posed; 3D /
  depth-resolved single-angle inversion is fragile and unnecessary here. One reflection measures one
  scalar projection **G·u** of the displacement field; twin selectivity makes the **amplitude a
  twin-domain map** ((202) lights one variant, (131) the other).
- **Photons are the currency and beam damage caps the dose** — not source brightness. NbN absorbs only
  ~1.5×10⁻⁴ of 10 keV photons, so damage is event-cascade (N-vacancy) driven; the empirical
  single-region ceiling is ~1h50min, and David's overnight macro already revisits a fresh region every
  ~39 min (≈2× safety margin). Tails run at tens of photons/pixel. This makes Poisson-aware /
  prior-regularized methods the correct family, and information-per-photon the figure of merit.
- **Bragg geometry breaks the simple transmission picture:** the probe **footprint** elongates 1/sinθ
  along scan-x (**5.15× at θ=11.2° for (202)**, 3.14× at θ=18.6° for (131)) — so with a ~150 nm spot
  (David's email) the *map resolution* is ~150 nm in y but ~0.5–0.8 µm in x. The sample-frame *pitch*,
  by contrast, is **isotropic ~20 nm — the sinθ pre-scaling of zpx is this experiment's universal
  convention**: every logged fly2d uses `zpx, ±X*sin(radians(samth.VAL))`, including **scan 289 on
  this very sample at the (131) angle** (zpx ±1.914 = ±6·sin 18.6°, notebook screenshot), and the
  data agree (speckle decorrelates **1.7–2.0× slower** along fast/x — measured, §3 — as an elongated
  footprint predicts and an elongated pitch does not; under the "5× coarser x-pitch" reading the two
  axes would decorrelate *equally*). One-line confirm from David that `Night5Ptycho()` does the same
  closes it.
- **Two different "elongations", easy to conflate — keep them apart.** The *surface footprint* stretches
  as 1/sin(incidence), and that is what limits **map resolution** (≈150 nm × 0.5–0.8 µm). The *far
  field*, however, sees that footprint projected onto the plane normal to the **exit** beam, so the
  aperture that sets the speckle scales as sin(exit)/sin(incidence) — and these are **asymmetric**
  reflections: 2θ = 47.1° with 11.2° in / 35.9° out for (202) → factor **3.02**; 55.9° with 18.6°/37.3°
  for (131) → **1.90**. Neither equals the 5.15 / 3.14 of 1/sin(incidence) alone. The detector samples a
  tilted Ewald slice; the effective probe must be "distilled" from the
  geometry.
- **The object confounds three things** — twin amplitude, strain phase-ramps (a strain gradient moves
  the peak inside the detector ROI → "phase-ramp pileup"), and mosaic tilt (the peak walks in χ as θ
  rocks). Separating them is the core inversion subtlety; a single fixed-θ COM gives only 2
  numbers/pixel and **cannot cleanly split tilt (2 DOF) from strain (1 DOF)** — amplitude *is*
  separable via intensity.
- **Structural limit — narrower than first stated:** the two *shared* 20 nm ptycho scans (309/313) are
  on different regions, so their cross-peak per-pixel fusion is impossible — statistical comparison
  only. But the scan sheet shows the campaign holds same-region cross-channel data we don't have yet:
  the (131) ptycho scans 312–315 share their regions one-for-one with **specular (111) ptycho scans
  316–319** (same zpy offsets, same night), and **rocking-map stacks 291–304 cover BOTH peaks ×7
  angles on a single region** of this film. High-res two-*twin* fusion remains a next-beamtime design
  item (fiducial-registered, ≤2 peaks per region's dose budget); twin+specular per-pixel work and
  multi-angle separation need only a Globus pull.

## 3. Results on the real data

Ran a first-look pipeline on both shared files. Key results and data-reduction facts:

- **Detector bad-pixel sentinel = 2³²−1** (3.58% of pixels incl. module-gap rows). Must be masked or
  every statistic is ~8 orders of magnitude wrong. `hdf5plugin` required (bitshuffle).
- **Raster grid = 71 × 41** (zpy slow × zpx fast; not 41×71 — resolved by map coherence and by
  COM-smoothness). Nominal 20 nm steps → ~1.4 × 0.8 µm field if the pitch is isotropic (see §2).
  **No serpentine artifact.**
- **No per-frame positions in the shared files** (`/entry/scan` is empty — David's email describes
  `Positioner/ZPX,ZPY` arrays that these exports don't contain). This is now a **soft** blocker for
  reconstruction, not a hard one: with the sinθ convention confirmed (§2), the scan center (CSV: zpx
  −0.239, zpy stepped per scan), 71×41 points, and 20 nm pitch, a **nominal position grid can be
  built from the scan sheet**, and joint position refinement (Calvo-Almazán 2019) absorbs the residual. The
  exact ranges sit in the beamline **Diffraction User Data Explorer** table (per-scan M1/M2
  MIN/MAX/NP — visible in David's notebook screenshots) → ask him to screenshot rows 308–313, easier
  than a master file. Still need **detector pixel size and beam center** (the explorer names the
  detector `s26_eiger2_cnm` → Eiger2, 75 µm px — confirm beam center) and the distance (scan sheet
  logs R=522.5; notebook elsewhere says 395/400/500 mm — which applied to 308–313?). Exposure was
  0.5 s/frame (notebook: "202 peak at 0.5 second exposure").
- **Energy ≈ 10 keV (9.6–10.5), and definitely NOT the 12.5 keV in the metadata.** The scan sheet logs
  the diffractometer angles (scan 309: samth 11.2°, twotheta 42.2°, gamma 23.15°; scan 313: samth 18.6°,
  twotheta 52.4°). With cos2θ = cos(tth)·cos(γ) both peaks index to standard NbN (a=4.39 Å, N=8 and
  N=11) at 10.0 keV. **Caveat found on review — an earlier "exactly 10.00 keV" was over-precise:**
  `twotheta`/`gamma` are detector-*arm* positions, but the peaks land 2.3–2.7° off the detector centre
  (pixels (813,591) and (854,577) vs a centre near (531,514)), so the true scattering angles differ from
  the logged ones by that much and the energy is only pinned to ~±0.5 keV until the beam centre is known
  (ask 3). Unaffected: 12.5 keV would need 2θ ≈ 37° against a logged 47°, a 10° gap far outside the
  ±2.7° the detector subtends — so it is excluded, the metadata field is stale, and the (202)/(131)
  agreement survives because both peaks are offset in the same sense.
- **Overnight dose management verified from the scan sheet:** scans 308→313 start ~39 min apart with
  zpy stepped 10 µm per scan — fresh region every scan, ~2× under the ~1 h 50 min damage ceiling.

**The physics that came out (see `twin_comparison.png`; all numbers reproduce with default
`quicklook.py` flags, ROI = peak ±64 px):**
- **(202)** is bright (median 229k counts/frame, ~14 ph/px in the ROI) with nearly flat integrated
  intensity (std/med 3.0%) but rich, coherent speckle — matching David's description. Peak-COM maps
  show **coherent lattice-tilt/strain domains** of hundreds of nm with boundaries sharp along the slow
  axis (along x, features are footprint-smeared — §2). Raw adjacent-frame speckle correlation
  0.90 fast / 0.89 slow.
- **(131)** is ~4.2× dimmer (55k counts/frame, ~3.4 ph/px) with ~2× stronger intensity modulation
  (std/med 6.5%) organized in diagonal few-hundred-nm bands, plus its own tilt/strain domains. Raw
  speckle correlation 0.66/0.65. **Attribution:** the stripe morphology itself is David's prior
  observation — his fly maps of this film (scan 289 "is this real?? yes it's real", 305–307, and the
  zoomed variants) already show the stripe domains at larger FOV. Our 20 nm step scan cross-validates that
  morphology; the new contributions are the **COM tilt/strain domain maps**, the **quantified
  modulation statistics**, and the corrected-correlation result below.
- **The correlation contrast is *pure photon statistics*.**
  Correcting the Pearson normalization for shot noise (Poisson noise variance = per-frame mean, so
  r_corr = cov/√((var−mean)(var−mean))) sends **both** reflections to the *same* underlying structural
  correlation: (202) 0.998 fast / 0.987 slow, (131) 0.999 / 0.983. The raw values are reproduced to
  ~1% by r_true × signal-variance-fraction (0.90 and 0.66 predicted from 14 vs 3.4 ph/px). Same
  structural continuity in both twins; the *information per frame* differs, which is the
  information-per-photon framing behind Plan B. (The shot-noise correction itself is standard
  XPCS/speckle-visibility practice; what is new is using it to separate structure from photon
  starvation across twin reflections.)

**Speckle geometry — the probe and the acquisition, measured from the data**
(`speckle_geometry.py` → `speckle_geometry.png`; inversion validated on synthetic speckle to <1%):
the speckle autocorrelation is the Fourier fingerprint of the coherently illuminated volume, so it
measures lengths David is still waiting on from the beamline.
- **The detector spends its pixels very unevenly — an acquisition-design result.** Speckle is 8.5 px
  FWHM on one axis but **1.4–1.7 px on the other, i.e. below the 2 px Nyquist floor**. In ptychographic
  terms the per-pattern field of view is λR/px = 864 nm, giving oversampling of **13.5× on one axis and
  only 1.6× ((202)) / 2.2× ((131)) on the other**. Consequences: rung 2 needs a
  probe-support constraint on the tight axis; the oversampled axis can be **binned 4–6× without
  loss** (cutting the ~300 GB data volume); and next beamtime should move the detector back or prefer more
  symmetric reflections. This is the Fisher/information-per-photon theme applied to the detector
  distance rather than the rocking angle.
- **The illumination has two length scales, and the probe size is confirmed at ~158 nm.** The
  autocorrelation is *not* a single Gaussian — fitting one is the wrong model and returns ~63 nm, an
  artefact lying between the two real scales (a two-component fit cuts the residual by 29% / 50%).
  The correct decomposition, agreeing across both scans:
  a **sharp** ACF peak → **157 nm ((202)) / 159 nm ((131)) illumination envelope**, and a **broad**
  pedestal → **~50 nm fine structure inside that envelope** (weights ≈ 42/58 and 46/54).
  The envelope number is confirmed three independent ways: the ACF sharp component, David's ~150 nm
  metadata spot size, and the real-space speckle decorrelation length measured from probe translation
  (161 nm / 146 nm). That the sharp component is anisotropic (σ 1.5 × 0.55 px) rules out detector
  charge sharing, and its width beyond the single zero-lag pixel rules out shot noise.
- **Consequence: the probe is structured, so rung 2 needs multi-mode probes.** A ~50 nm structure
  inside a ~158 nm envelope is what a zone plate with a central stop produces (ring structure), and is
  equivalently what a ~50 nm transverse coherence patch inside a ~158 nm beam looks like. Either way
  the correct algorithmic response is **mixed-state / multi-mode ptychography** (Thibault & Menzel
  2013) rather than a single coherent probe.
- **The discriminating experiment is already collected.** If the 50 nm scale is sample-driven (a
  coherently diffracting mosaic block, which would make it a *science* observable), it must change on
  the Au standard with its large EBSD-known domains, and across the 6 nm / 40 nm films with their very
  different domain scales. If it is probe-driven it will not move. So the Au data serve as a
  physics test as well as solver bring-up, and need only a Globus pull.
- ~~Earlier claim, now retracted:~~ an initial single-Gaussian fit reported "the coherently diffracting
  region is 2–5× smaller than the illuminated footprint." That gap was entirely an artefact of the
  wrong model; with two components one of them *is* the footprint. `speckle_geometry.py --selftest`
  now reproduces this exact failure mode on synthetic data and asserts against it.

This covers the group's "step 0" asks — intensity maps, COM maps, adjacent-pattern correlation,
anisotropy — on real data with a tested pipeline, plus the noise-corrected correlation.

## 4. Directions (one staged ladder, not competing tracks)

Everything lives on one pipeline, and each rung de-risks the next:

0. **Reconstruction-free maps** — intensity, COM, virtual-aperture (satellites), speckle-correlation
   (raw **and shot-noise-corrected**) + decay-length **anisotropy**. *Done for (202)/(131).*
   Pedigree: Chahine K-map/XSOCS, Kriegner
   (twin mapping in Bi₂Te₃ — same ABC/ACB degeneracy, EBSD-validated), Morgan **Ptychographic X-ray
   Speckle Tracking** (the exact precedent for "correlate adjacent patterns → map"; reconstruction-free,
   open code, and **no Bragg implementation exists yet** — a methods gap).
1. **COM → tilt/strain decomposition** — *Done (`gu_proxy.py` → `gu_proxy.png`); the data revise the
   original premise.* The plan was Shack-Hartmann-style integration to a smooth G·u map;
   the measurement says that shortcut mostly doesn't exist here: Helmholtz-decomposing the domain-scale
   COM field at the most-favorable detector→sample frame, **at most 44% ((202)) / 34% ((131)) of the
   field power is compatible with *any* scalar strain-gradient surface — the majority is rotational**,
   i.e. the mosaic is twist-dominated in real space (the nanoscale counterpart of David's ensemble
   "twist mosaicity"; robust across discrete-pairing, rotation-fit, and smoothed variants). Deliverables:
   the **twist-domain map Ψ** (new modality, hundreds-of-nm twist texture; weakly anticorrelated with
   intensity, r≈−0.24/−0.25 in *both* scans — suggestive, preliminary) and the **gradient-part G·u proxy**
   (±25 pm scale, provisional calibration) which still serves as solver initialization for its share of
   the field. Consequence, stated quantitatively: single-angle COM cannot be summarized as strain —
   the 3-angle separation (Holt PRL 2014) and full phasing (rung 2) are *necessary*, not optional.
   **And the separation is executable on existing data:** the scan sheet shows same-region rocking-map
   stacks through both peaks for *all three thicknesses* — 10 nm: scans 291–304 (7 angles × 2 peaks,
   0.3° steps); 6 nm: 326–363 (~10 angles × 2 peaks, two detector distances); 40 nm: 369–384
   (8 angles × 2 peaks). A thickness-series twist/strain separation needs a Globus pull, not beamtime.
2. **Single-peak 2D BPP reconstruction** — now with three measured design constraints from §3:
   a **probe-support constraint** on the undersampled detector axis (oversampling only 1.6–2.4×),
   **mixed-state/multi-mode probes** (the ~50 nm structure measured inside the ~158 nm envelope), and
   the exit-projection aperture factor (3.02 / 1.90, not 5.15 / 3.14) in the geometry layer.
   Original plan: geometry-distillation preprocessing + effective elongated
   probe + a 2D engine. No off-the-shelf Bragg engine exists (pty-chi has zero Bragg support; EASYBRAGG
   is consortium-gated 3D), so the shortest path is the small, well-documented geometry layer in front of
   pty-chi/LSQML. Precedents: Hruszkewycz 2012/2014, and **Burdet 2017 did BPP on a Nb film** — nearest
   published system; expect a refinement stack (positions, ramps, peak centering).
3. **PtychoPINN-Bragg and information-optimal acquisition** (see §5, Plan B).
4. **Multi-peak u(r) fusion** — next beamtime only (Kandel 2023 multi-peak BPP; needs same-region
   fiducial registration).

## 5. Plans

### Plan A — Reconstruction-free twin & strain mapping on the shared data
- **Deliverables exist:** `twin_comparison.png` (2×4: intensity, COM, raw corr, **noise-corrected
  corr**), `gu_proxy.png` (rung-1 Helmholtz/twist decomposition), `speckle_geometry.png`
  (probe/coherence + detector sampling), `quicklook.py`, `speckle_geometry.py` (both with selftests), `data/*_maps.npz`.
- **Summary:** the PI's "4 numbers per image" idea on real data, extended twice: (1) the raw correlation contrast between twins is fully explained by photon statistics
  (corrected r identical at ~0.98–0.99), so the correlation map is a photon-budget-aware structural
  probe; (2) the COM field, decomposed, is **twist-dominated** (≤44%/34% strain-gradient-compatible)
  — a real-space twist-domain map that echoes David's ensemble twist mosaicity and *quantifies why*
  phasing + 3-angle scans are necessary. It also identifies a methods gap (Bragg speckle tracking).
- **Optional extensions:** quick NMF domain segmentation; per-frame speckle-orientation
  (structure-tensor) map — the peak's internal streak texture (visible in David's notebook detector
  images) carries local fault-plane orientation, and this is the PI's original "rotate so the lines
  are parallel" idea, which no beamline tool provides. Virtual-aperture maps are not new: David's Data
  Explorer GUI already makes ROI-box darkfield/satellite maps live at the beamline (his satellite
  maps came from it), so a scripted version adds reproducibility, not capability.
- **Caveats:** grid orientation inferred (71×41, validated by coherence + no-serpentine checks);
  pitch convention pending David (likely isotropic 20 nm — macro evidence); COM = mixed tilt+strain,
  not clean strain; map resolution is probe-limited ~150 nm and 1/sinθ-elongated along x;
  (202)/(131) are different regions → statistical comparison only.

### Plan B — PtychoPINN-Bragg + information-optimal acquisition
- **B1 — amortized Bragg inversion.** PtychoPINN already has what every existing DL-Bragg effort
  (AutoPhaseNN, NanobeamNN, DONUT) lacks: the **overlap constraint + Poisson loss**. Repo audit says the
  forced changes are modest (effective anisotropic probe via existing custom-probe ingestion; `nongrid`
  coords already supported; complex object already; `ci_intensity_v2` Poisson profile). **Seed
  figure:** `photon_budget.png` — radial counts/pixel of the real scans (10 ph/px by r≈50 px, 1 ph/px by
  r≈105 px; (131) an octave lower) — the count target any forward simulation / Poisson loss must hit.
  The `quicklook.py` selftest already contains the twin-domain forward model to extend.
- **B2 — Fisher-information acquisition design.** Answers meeting-1's literal "where on the rocking
  curve?": for a Gaussian rocking curve, information about strain is **zero at the peak and maximized at
  Δ\* = σ√2 ≈ 1.41σ off-peak**. **Seed figure:** `fisher_rocking.png` (with the σ√2 optimum
  marked; caption holds step-size/overlap fixed). The noise-corrected correlation in §3 shows
  photon-statistics accounting changing conclusions on this dataset, and the **detector-sampling
  result** (§3, oversampling 13.5× vs 1.6×) is the same argument applied to a second acquisition
  knob, with a direct recommendation. Together they make the
  information-allocation case concrete: rocking offset, detector distance, and pixel binning are all
  choices about where the photons' information goes.
- **Relation to ANL's ML work (Cherukara/Zhou/Holt — DONUT/NanobeamNN):** complementary, and already
  routed through 26-ID. Their methods are per-pattern geometric regression; PtychoPINN adds phase
  retrieval with the overlap constraint and a Poisson loss.

**Recommendation:** A first (real data, the group's stated ask), B in parallel (needs no new data).
Reconstruction (§4 rung 2) follows and is no longer hard-blocked: nominal positions can be built from
the scan sheet (§3), David's Data Explorer rows make them exact, and joint position refinement
absorbs the rest.

## 6. Questions for David
1. Per-frame **positions**: master file, **or just a screenshot of the Diffraction User Data Explorer
   rows 308–313** (the M1/M2 MIN/MAX/NP columns) — that plus the confirmed sinθ convention gives us
   exact position grids without any file transfer. (The shared exports have an empty `/entry/scan`.)
2. One-line confirm: **`Night5Ptycho()` pre-scales zpx by sin(samth)** like every fly2d in the sheet,
   and zpy is the 71-point axis. (Exposure 0.5 s already confirmed from the notebook.)
3. Detector **beam center + distance for scans 308–313** (detector is `s26_eiger2_cnm` → Eiger2,
   75 µm px). The R column decodes as nominal+~22.8 mm (423.1↔"400 mm", 522.5↔"500 mm"), so 308–319
   ran at the **nominal 500 mm** setting — confirm whether the true sample–detector distance is 500 or
   the 522.5 readback (4.5% on all q/pm calibrations).
4. **Energy metadata is stale:** logged angles index both peaks near 10 keV (±0.5 keV until the beam
   centre is known, §3); `DCM Energy` reads 12.5 keV and should be corrected in the master file.
4b. **Probe structure:** the speckle autocorrelation resolves a 157/159 nm envelope (matching the
   ~150 nm metadata spot and the measured decorrelation length) *plus* ~50 nm structure inside it (§3). Useful to know: which zone plate and
   defocus ZP2 was at for 308–319, whether it has a central stop / order-sorting aperture (which would
   produce exactly this ring structure), and the beamline's transverse-coherence estimate (slit and
   source sizes). That settles probe-structure versus partial coherence, and either way tells us
   how many probe modes rung 2 should solve for.
5. Status of 308/312 + a Globus path to watch.
6. **One gold-standard file for rung-2 solver bring-up** (per meeting 2's "start with the Au
   standard" — EBSD-known domains, e.g. Au 020 step scans 247–249 — and meeting 1's "start from the
   111 with the working recon": a specular 111 ptycho scan, e.g. 102/104, serves the same role).
   The NbN 10 nm files suit the reconstruction-free analysis; the *reconstruction* step needs
   validation targets with known structure first.

**Globus pull priorities** (from the full scan-sheet inventory; sizes are per-scan ~0.7 GB step /
smaller fly):
1. Au standard: 247–249 (+ Au rocking maps 236–243) — solver bring-up.
2. **291–304** — same-region **two-peak × 7-angle** rocking maps, 10 nm film → Holt tilt/strain
   separation + first cross-peak registration.
3. **316–319** — specular (111) ptycho, *sharing regions with (131) scans 312–315* → per-pixel
   twin-selective vs twin-insensitive two-channel analysis; 314/315 for (131) statistics.
4. 326–363 (6 nm) and 369–384 (40 nm) two-peak rocking stacks → thickness-series twist/strain story.
5. 179–189 — 6 nm specular ptycho (ANL's "start from the 111" dataset).

## 7. Key references (full list in analysis notes)
- **Methods (Hruszkewycz arc):** Nano Lett 2012 (first nanofocused BPP) → **PRL 110, 177601 (2013):
  BPP phase → picometre displacement/polarization map of 180° nanodomains in an epitaxial film — the
  direct template for our twin-domain displacement deliverable** → Holt PRL 2014 (tilt/strain
  separation from BPP phase) → Nat Mater 2017 (single-angle 3DBPP — the fragility argument for staying
  2D) → PRA 2016 (BPP domain sensitivity) → Kandel Opt Express 32, 19594 (2024) + Porter JAC 2026
  (multi-peak BPP → vector u(r), robust inter-peak coupling) → Maddali JAC 2020 (shear correction).
  **Burdet PRB 96, 014109 (2017) — BPP on a Nb film, Robinson group (NOT Hruszkewycz): nearest
  experimental precedent, correctly attributed.**
- **Experiment design / registration (his corpus answers our two blockers):** Calvo-Almazán Sci Rep 9,
  6386 (2019) — jointly reconstruct the crystal AND each frame's true angle, i.e. the software fix for
  our missing/uncertain scan positions; Zhou PRB 2024 (Takagi–Taupin for focused-nanobeam Bragg — flags
  when kinematic BPP phase is biased by the dynamically-diffracting substrate beneath an epitaxial film).
- **Reconstruction-free:** Chahine JAC 2014 (K-map) + XSOCS, **Kriegner JAC 2017 (twin mapping, our exact
  degeneracy)**, **Morgan JAC 2020 (X-ray speckle tracking)**, py4DSTEM, Venderley PNAS 2022 (X-TEC).
- **Dose / information:** Wei PRA 2020 (Poisson CRB for ptychography), Bouchet Nat Phys 2021
  (max-information states), Dwyer-Paganin PRB 2024 (4D-STEM Fisher info), Hegerl-Hoppe dose fractionation.
- **ML-CDI:** PtychoNN (APL 2020), AutoPhaseNN (npj 2022), NanobeamNN/DONUT (2025, the 26-ID ML frontier —
  stops short of phase retrieval), PtychoPINN (Sci Rep 2023).
- **Materials:** Chang SuST 2021 (dense rotational twins in δ-NbN(111)), Wright PRMat 2021 (growth-mode
  control), Yan Nature 2018 + Kim Commun Mater 2021 (why: qubits/SNSPDs).
- **Novelty:** Hruszkewycz's corpus contains every ingredient — nanofocused BPP on
  thin films, phase→displacement, defect/stacking-fault phasing, angular-registration recovery,
  multi-peak vector imaging — but **no NbN or superconducting-film paper**; this would be the program's
  first application to superconducting quantum-device materials.
