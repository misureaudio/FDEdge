# Sliding-window sweep for fdedge on long signals (audio)

**Status.** Written 2026-10-17 from the v1.0.0 discussion
("how should the sliding window be implemented to drive fdedge local
detection of transient-related edges in a large wav/flac file?").
Every number below was measured headless on base MATLAB R2024a against
the v1.0.0 library. Reference implementation: `examples/sweep_demo.m`
(self-contained, synthesizes its own 60 s test file and scores itself).

Companion to `docs/real_signal_usage.md` (single known-edge window);
this document covers the **unknown-edge** case: finding where the
edges are in a long file.

## 1. Design in one paragraph

Sweep the file with a window of length **M** (the length the detector
is calibrated for, 1024 samples = 23.2 ms @ 44.1 kHz in the verified
runs) at hop **H = M/4** (75% overlap). Per window: subtract the
pre-edge baseline (median of the left quarter), smooth, one
Caputo–GL evaluation at a fixed detection order (α = 0.70). The
**detection statistic is the step-height profile** — the boxcar
level difference of the GL response across each trial edge position
(computed in O(M) via prefix sums) — *not* the slope. Threshold the
per-window max, NMS in window units, dedup in sample units. Each
surviving candidate then runs the **full v1 pipeline**
(`Tracker.update` with the α-bank + sub-sample refine) and must pass
the two-tier gate (stepSNR > 3, position sanity).

## 2. Why per-window GL (and not whole-file)

The Caputo–GL response at a point depends on the **entire window**
(non-local fractional kernel), so a whole-file evaluation is both
O(N²)-class and physically wrong: the response at sample *i* would be
contaminated by everything before *i* (the finite-window boundary
artifacts of pitfalls P3/P5 apply at the file edges, and slow drift
acts like the DC-tail of §3.1 of the companion doc, on a global
scale). Per-window evaluation makes the operator **local with memory
M**: O((N/H)·M²) total, O(N) memory (a 24 h mono 44.1 kHz file is
~2.5 GB as doubles — fits in RAM). The per-window median baseline
removes DC and slow drift locally.

## 3. Stage 1 — coarse scan (one GL per window)

Per window (M = 1024, H = 256, α = 0.70, s = 8):

```matlab
w  = y(n0:n0+M);
w  = w - median(w(1:round(0.25*(M+1))));   % pre-edge baseline
fs = smoothGauss(w, 8);
R  = caputoGL(x, fs, 0.70);
% step-height profile (boxes of Wbox = 96 samples, gap dgap = 32):
cR = [0; cumsum(R)];
ii = (dgap+Wbox) : (M+2-dgap-Wbox);
Hv = (cR(ii+dgap+Wbox) - cR(ii+dgap) - cR(ii-dgap) + cR(ii-dgap-Wbox+1)) / Wbox;
```

`H(p)` = mean response in `[p+32, p+128]` minus mean in
`[p−128, p−32]`, restricted to the interior zone `[0.25M, 0.75M]`
(the left quarter is baseline context, the boxes must stay interior).

### 3.1 Why step height, not slope (measured)

The slope (12-point central difference of R, the PoC's *localization*
metric) is the **wrong detection statistic**: its noise floor sits
above the edge signature of weak onsets.

| statistic | noise (23 windows over 60 s, per-window max) | weakest onset (A = 0.3, β = 0.55) | margin |
|---|---|---|---|
| max\|slope\| | median ≈ 25, q99.9 ≈ 44 | ≈ 25 | ~1× (undetectable) |
| max step height H | median 0.072, q99 0.126 | ≈ 0.37 | ~11σ |

Step height works because the GL of a one-sided cusp at an order
near its exponent is a **step of height ≈ A** in the response; the
boxcar difference averages 96 samples of each side (σ of the
difference ≈ 0.033 at σ_noise = 0.02), while the slope averages only
12 and amplifies noise by M/12.

### 3.2 Localization: centroid of the step profile (measured)

For a true edge the profile is a **trapezoid symmetric about the
edge** (rise over [T−128, T−32], flat plateau [T−32, T+32], fall over
[T+32, T+128] in sample units), so its positive-part centroid
(windowed to ±154 samples around the argmax) equals the edge within
a few samples:

```matlab
Hc = max(Hv(lo:hi), 0);  xc = ii(lo:hi); xc = xc(:);
xloc = (sum(xc .* Hc) / max(sum(Hc), eps)) / M;   % note the /M!
```

The slope argmax is a **bad locator** for β > α edges: their
signature in the response is a rising ramp (R ∝ d^(β−α) from 0), whose
edge slope is only ~1σ above the slope noise floor for moderate A —
in the verified runs it missed the β = 0.85, A = 0.5 onset by 80
samples (0.078 in window units), beyond the refine's ±0.02 search
window. The centroid made the same onset land within 14 samples.

### 3.3 Threshold, NMS, dedup (measured)

- `thrH = 0.2` (≈ 1.6× the q99 noise max 0.126). Weakest verified
  onset: 0.37.
- NMS in **window units** (min gap `M/(2H)` = 2 windows), then a
  second dedup in **sample units** (keep the best `hmax` within
  `0.3·M` = 307 samples) — one edge sits in the interior zone of up
  to two adjacent windows and must not be reported twice.
- Measured on the 60 s test file: 24 raw windows → 8 after window
  NMS → 8 after position dedup (the 8 onsets).

## 4. Stage 2 — per candidate (the full v1 pipeline + gate)

```matlab
cfg   = config('Profile','beta','M',M);
calib = gridSlopePeak(cfg);
trk   = Tracker(cfg, calib);
[st, ~] = trk.update(w, xguess);          % xguess from the centroid
Th = (n0-1) + st.xraw*M;
R2 = caputoGL(x, smoothGauss(w, cfg.sGrid), st.braw);
iE = min(max(round(st.xraw*M), 1), M);
stepSNR = mean(R2(iE+32:min(M+1, iE+128))) / std(R2(1:64));
Ahat    = mean(R2(iE+52:min(M+1, iE+205)));
ok = (stepSNR > 3) && (abs(st.xraw - xguess) <= 0.0201);
```

- **stepSNR > 3 is the real gate.** In the verified runs every one of
  the 23 noise-driven candidates scored ≤ 1.4 (the step at a noise
  location is a difference of two noise means ≈ 0); the weakest true
  onset scored 7.8. The margin is 5×+ in both directions.
- **Position sanity ≤ 0.0201** is the refine's own search width (the
  refine returns a position inside `xguess ± 0.02` **by
  construction**, so this never rejects a refine result — it is a
  diagnostic for pipeline misuse, *not* a false-alarm gate). A
  locator failure (xguess 0.05–0.08 off) shows up as the refine
  pinning to the search boundary **plus** a low stepSNR at that
  location — the stepSNR gate is what rejects it (measured: the
  44.2 s onset candidate pinned 80 samples off scored 2.73 < 3).
- `Ahat` (response level 52–205 samples past the edge, where the
  response sits at ≈ A·Γ-ratio for α ≈ β) is the event amplitude
  estimate.

## 5. Verified performance

Test file: 60 s @ 44.1 kHz, 8 one-sided power-law onsets
(A = 0.9, 0.6, 0.4, 0.8, 0.3, 1.0, 0.5, 0.7; β = 0.55, 0.70, 0.70,
0.85, 0.55, 0.70, 0.85, 0.70; M-sample attack + 100 ms exponential
decay) + DC 0.3 + 50 Hz hum (0.01) + 0.1 Hz drift (0.05) + white
noise σ = 0.02.

**Main run (rng(99)): 8/8 detected, 0 false alarms.**

| k | dT0 (smp) | dTref (smp) | β true | β̂ | A true | Ā | stepSNR |
|---|---|---|---|---|---|---|---|
| 1 | −16 | +1 | 0.55 | 0.5622 | 0.9 | 0.841 | 24.5 |
| 2 | 0 | −13 | 0.70 | 0.7211 | 0.6 | 0.618 | 9.5 |
| 3 | −12 | −31 | 0.70 | 0.5500* | 0.4 | 0.301 | 13.3 |
| 4 | +14 | +3 | 0.85 | 0.5500* | 0.8 | 0.463 | 8.7 |
| 5 | −12 | −18 | 0.55 | 0.6079 | 0.3 | 0.318 | 13.2 |
| 6 | −3 | −3 | 0.70 | 0.6607 | 1.0 | 0.888 | 21.4 |
| 7 | +6 | +5 | 0.85 | 0.7162 | 0.5 | 0.405 | 7.8 |
| 8 | −1 | +1 | 0.70 | 0.6806 | 0.7 | 0.651 | 24.4 |

dT0 = centroid error before the refine; dTref = final position error
(≤ 31 samples = 0.7 ms). \* β̂ pinned at the bank minimum for
β = 0.85 edges at moderate A — the v1 slope-peak β metric is
weakest for β ≫ α (ramp signature); position and step gates are
unaffected. Known v1 limitation, not a sweep defect.

**5-seed robustness (rng(1001..1005)): 39/40 detected, 0 false
alarms in all 5 seeds, max |dTref| = 47 samples (1.07 ms).** The
single miss was the A = 0.4, β = 0.70 onset (weakest class) on one
seed — a gate-boundary case (stepSNR just under 3).

**Cost:** stage 1 = 1.4 s, stage 2 = 0.3 s, total 1.7 s for 60 s of
audio ≈ **35× realtime** on this machine (i7-class, R2024a base).
A 1 h file ≈ 1 min. Memory O(N) for the file + O(M) per window.

## 6. Model caveat — the cusp is a local *attack* model

Discovered while building the test file: a synthetic one-sided power
law that **keeps rising forever** makes the detector flag the
ongoing rise as (several) edges — correctly, because the response
there is a real power-law response. Real audio transients are
**bounded**: attack (power law over ~10–30 ms) followed by decay or
sustain. The sweep finds the **attack corner** (the onset) — which is
exactly what transient-edge detection wants — and ignores the decay
tail (a smooth exponential has no corner: its GL response is flat-to-
ramp, and the step profile of a flat response is ≈ 0).

If your signal contains genuinely sustained level shifts (e.g. a
permanent DC offset change), the detector will report the shift edge
once — as intended.

## 7. Tuning knobs

| knob | verified value | effect of raising |
|---|---|---|
| `thrH` | 0.2 | fewer candidates, weaker onsets missed (weakest verified 0.37) |
| `gateSNR` | 3 | fewer false alarms (noise max measured 1.4), weaker onsets rejected (weakest verified 7.8) |
| `M` | 1024 (23.2 ms) | longer attack coverage; rebuild `calib` per M (D7); P4 guard M ≤ 8192 |
| `H` | M/4 | smaller = denser coverage, linear cost; the [0.25, 0.75] interior zone + dedup make M/4 sufficient |
| `Wbox`, `dgap` | 96, 32 | box size vs. position resolution of the step profile; 96/32 gives a 64-sample plateau (±32) around the edge |
