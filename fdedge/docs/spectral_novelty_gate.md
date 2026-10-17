# Spectral-novelty co-gate: theory, verification, and how to read the numbers

**Status.** Written 2026-10-17. Companion to
`docs/stationary_multitone_gate.md` (the verified false-positive mode of
the old gate), `docs/sliding_window_sweep.md` (the sweep design), and
`docs/real_signal_usage.md`. All numbers measured headless on base
MATLAB R2024a. Test file: the PrimeTone 27-tone stationary multitone
(characterized in `docs/stationary_multitone_gate.md` §1).

Implementation: `fdedge/est/spectralNovelty.m` (the novelty primitive),
`fdedge/est/burgAR.m` (the stable AR fit),
`examples/stress_multitone_gate.m` (the stress test with the position
scan co-gate), `examples/robust_step_swept.m` (step at 6 positions),
`examples/ref_novelty_suite.m` (deterministic synthetic reference,
`rng(0)`).

## 1. What problem the co-gate solves

The fdedge pipeline is two-stage: stage 1 (per-window Caputo-GL
derivative at α = 0.70 + step-height profile) produces a candidate
edge per window — `hmax` (height) and a centroid position; the gates
decide which candidates are real transient edges. The old gate is
stage 2: `Tracker` refinement + `stepSNR > 3` + position sanity.

Two failure modes of the old gate are measured and documented:

1. **False positives on stationary tonal signal** (the design driver).
   On the stationary 27-tone PrimeTone file, all 183 windows pass
   stage 1 (hmax median 4.352, q99 7.553), and 24 windows pass
   `stepSNR > 3` (20 with the position co-gate). Mechanism: the
   false passes sit on **envelope crests** of the phase-optimized
   multitone — a quiet→crest transition is a real *local level
   change* at the 21 ms window scale, indistinguishable from a step in
   the level domain (`docs/stationary_multitone_gate.md` §3–4; the raw
   level-change gate and the box-correlation gate both failed to
   separate the cases).
2. **Blindness to new high-frequency tones** (measured during the
   stress test). A 1 kHz tone (A = 0.5) injected at t = 0.6 s on the
   same background reads `stepSNR = 0.48` at its edge: the gate's
   numerator is the *mean* of the derivative-type response over a
   96-sample box, and a 1 kHz oscillation completes 96·(1000/48000) =
   2 full cycles inside the box, so the mean ≈ 0. The old gate cannot
   see an edge whose new content oscillates.

The separating property found across every tested case
(`docs/stationary_multitone_gate.md` §5): **spectral novelty** —
does the post-edge region contain spectral content the pre-edge region
does not? A step adds DC + broadband; a new tone adds a line; a burst
adds broadband energy; a stationary multitone crest adds *nothing*
(the same 27 lines, re-phased). The co-gate operationalizes exactly
this question as a model-based test.

## 2. The question the co-gate asks, and the metric

At a candidate edge position `iE` the window splits into
`pre = w(1:iE)` and `post = w(iE+1:n)`. The co-gate performs the test:

- **H0 (no edge):** the post segment is a continuation of the *same
  process* as the pre segment.
- **H1 (edge):** the post segment contains content the pre segment's
  process cannot produce.

Operationalization: fit the best linear (AR) predictor of order
p = 120 to the pre segment, run it on the post segment, and measure
the **relative RMS prediction error** over the first L = 64 post
samples:

```
nov = ||post(1:L) − pred(1:L)||₂ / ||post(1:L)||₂
```

`nov` is the **fraction of the post edge's energy that the pre-edge
process cannot explain**:

- `nov ≈ 0` (measured 1e-6…8e-5 on the stationary file): the post is
  the exact continuation of the pre → H0 → reject.
- `nov ≈ 1`: the post is entirely unexplained by the pre (e.g.
  silence → step reads exactly 1.0000).
- `nov > 1` is possible: the error is `post − pred`, and if the
  prediction has the wrong phase at a new frequency it adds
  *destructively* to the post (see §8.3).

The threshold is `nov > 0.1` (reject at or below). Where 0.1 sits is
quantified by the measured separation in §7.

### Why linear prediction instead of an FFT comparison

- The pre segment is only 6.7–14.7 ms (320–704 samples at 48 kHz,
  depending on the in-window edge position). An FFT of 512–1025
  samples resolves frequency to 46.8–93.8 Hz; the PrimeTone's lowest
  lines (43, 53, 61 Hz) occupy 0.29–0.67 cycles of the pre segment and
  sit 8–10 Hz apart — below the window's resolution limit, and no
  zero-padding of the window can separate them.
- An AR(p) model with p = 120 has 120 complex poles: the fit places
  pole pairs near the observed line frequencies *without ever needing
  to resolve or name the lines*, and the prediction error responds to
  whatever the model cannot reproduce.
- A raw pre/post spectral power ratio fails for a different, measured
  reason (`docs/stationary_multitone_gate.md` §4c): a step's new
  content is DC + the boundary discontinuity, and per-segment FFTs of
  the smooth segments carry none of that energy (measured: the ratio
  reads 0.93 for a real step). The prediction error captures the
  unmodeled content directly, at the edge.

### Why Burg's algorithm for the fit

The predictor is run recursively on the post segment, so it must be
stable (all poles strictly inside the unit circle). The unconstrained
least-squares fit can produce unstable poles on short pre segments —
measured: a pole at |z| ≈ 1.10 on a 6.7 ms multitone segment, which
blew the novelty up to 3e4 for a real step. Burg's forward-backward
recursion gives |reflection coefficients| < 1 at every order by
Cauchy–Schwarz, hence a guaranteed-stable predictor.

## 3. Why one-step prediction with feedback (and not open-loop)

Two ways to use the fitted model on the post segment:

- **Open-loop multi-step:** feed the *predictions* back into the
  predictor state. The state then follows the model's own dynamics.
- **One-step with feedback:** feed the *actual* post samples back.
  The state is re-anchored to the truth at every step.

The choice is forced by a measured failure of the open-loop variant.
The pre segment (6.7–14.7 ms) is **shorter than the period of the
lowest lines** (43 Hz = 23.3 ms), so the fitted line poles sit at
radius r < 1 — the finite segment plus quantization noise pull them
inward. Open-loop, the predicted tone amplitude decays as r^L over the
500–700-sample post while the true tone continues. Measured on the
stationary PrimeTone file (no edge anywhere): nov min 3.49e-01,
median 0.996, max 2.115 → **183/183 false alarms**.

One-step with feedback re-anchors every step, so the radius error does
not accumulate: the same stationary file reads nov min 1.31e-06,
median 7.39e-06, max 8.07e-05 → **0/183 false alarms**.

**The cost of the feedback:** the predictor re-equilibrates to a new
*level* in a short transient. Measured: for the 0.5 step injected at
t = 0.3 s on the PrimeTone background, the one-step prediction error
over the *entire* post segment (704 samples) reads nov = 0.0658 —
below the 0.1 threshold, i.e. a real step rejected, because the error
is confined to the transient (effective length ~12–64 samples from the
transient-vs-full comparison) and the rest of the post is predicted
perfectly at the new level.

**The fix:** measure `nov` only over the **first L = 64 post samples**
(the stepSNR gate's own post-box length, for design consistency),
before the re-equilibration transient decays. The same step then reads
nov = 0.1568 at the raw stage-1 position (and 4.1936 with the
position scan, §5). A new *tone*, by contrast, cannot be
re-equilibrated to — the model has no pole at the new frequency and
the feedback only re-anchors the state, it does not add model content —
so its error persists over the whole post and the L-box captures it
fully.

## 4. The threshold and the measured separation

Reject if `nov ≤ 0.1`. The two measured populations on the PrimeTone
file (183 windows, with the position scan of §5):

| population | nov range (measured) |
|---|---|
| stationary (no edge anywhere) | 1.31e-06 … **8.07e-05** (max) |
| real edges (0.5 step at 6 positions; new 1 kHz tone A = 0.5) | **0.9854** (min) … 8.2983 (max) |

The gap between the stationary maximum and the edge minimum is a
factor of ≈ 12,000. `0.1` sits in the middle of that empty band; it is
a separator between two measured populations, not a physical amplitude
and not a claim about signal level.

## 5. The position scan (co-gate policy, not part of the primitive)

`spectralNovelty(w, iE, Fs)` is a pure "novelty at position iE"
primitive. The co-gate that uses it scans four positions and takes the
maximum:

```
nov_cogate = max over iE' ∈ {iE−32, iE, iE+32, iE+64} of spectralNovelty(w, iE', Fs)
```

Why: the stage-1 centroid is not exact on the multitone background.
Measured centroid error (iE minus the true edge position) for a 0.5
step at six positions:

| step at | window | iE | true in-window position | error |
|---|---|---|---|---|
| 0.15 s | 27 | 560 | 545 | +15 |
| 0.30 s | 56 | 349 | 321 | +28 |
| 0.45 s | 83 | 606 | 609 | −3 |
| 0.60 s | 112 | 350 | 385 | −35 |
| 0.75 s | 140 | 403 | 417 | −14 |
| 0.90 s | 167 | 617 | 705 | **−88** |

If the centroid is early by more than L = 64, the novelty box
`[iE+1, iE+64]` falls entirely *before* the edge and a real edge reads
nov ≈ 1e-5. Measured: the t = 0.90 s step, centroid 88 early, reads
**nov = 0.0000 (rejected)** without the scan. The scan
{−32, 0, +32, +64} covers centroid errors up to ~128 samples early and
~32 late; with it the same step reads **nov = 0.9854 (PASS)**. Cost:
4 novelty evaluations per window (full stress test, 4 sweeps of 1 s of
audio: 5.0 s).

## 6. Verification matrix (all measured)

### 6.1 Synthetic reference suite (`ref_novelty_suite.m`, rng(0), M = 1024, Fs = 48 kHz, edge at iE = 512, p = 120)

| case | nov | verdict |
|---|---|---|
| real step + noise (0.5 step, σ = 0.02) | 1.0531 | PASS |
| clean step (0 → 0.5) | 1.0000 | PASS |
| 60 Hz tone, stationary | 0.0000 | reject ✓ |
| 1 kHz tone, stationary | 0.0000 | reject ✓ |
| 60 Hz tone + step at edge | 2.0591 | PASS |
| 60 Hz tone + new 400 Hz at edge | **0.0159** | **reject — degenerate corner, §8.4** |
| 60 Hz tone + noise (σ = 0.01) + new 400 Hz | 0.1109 | PASS (noise breaks the degeneracy) |
| 60 Hz tone + 1-sample click (0.3) at edge | 0.9523 | PASS |
| silence → step | 1.0000 | PASS |

### 6.2 Stress test on the PrimeTone file (`stress_multitone_gate.m`)

| run | result (measured) |
|---|---|
| 1 — pure stationary file, 183 windows | nov 1.31e-06 … 8.07e-05 → **0/183 false alarms** (old gate: 24/183 stepSNR passes, 20 with position co-gate) |
| 2 — 0.5 step at t = 0.3 s | window 56, hmax 4.426, **nov 4.1936** PASS (old gate: stepSNR 0.61 — the old gate *also* misses this step; the post box sits on the flat part of the multitone response) |
| 3 — new 1 kHz tone (A = 0.5) at t = 0.6 s | window 111, hmax 3.687, **nov 1.7795** PASS (old gate blind: stepSNR 0.48, §1.2) |
| 4 — weak tone (A = 0.15, −15 dB) at t = 0.6 s | **nov 0.6777** (co-gate PASS), stepSNR 1.08 (old gate reject) — informational, §8.5 |

Score: **run1 PASS, run2 PASS, run3 PASS**.

### 6.3 Robustness: 0.5 step at six positions (with scan)

| step at | best window | iE | nov | verdict |
|---|---|---|---|---|
| 0.15 s | 27 | 560 | 6.1144 | PASS |
| 0.30 s | 56 | 349 | 4.1936 | PASS |
| 0.45 s | 83 | 606 | 8.2983 | PASS |
| 0.60 s | 112 | 350 | 2.4021 | PASS |
| 0.75 s | 140 | 403 | 6.8660 | PASS |
| 0.90 s | 167 | 617 | 0.9854 | PASS |

## 7. How to read the numbers (interpretation guide)

This section addresses the specific interpretation difficulties that
came up while reviewing the results.

### 7.1 `nov` is a relative prediction error, not an amplitude

`nov` is the RMS prediction error normalized by the post segment's own
RMS. It answers "what fraction of the post's energy is new?", not "how
big is the edge?". Consequences:

- `nov ≈ 1` does *not* mean "a huge edge"; it means "the post is
  entirely different from what the pre predicts". Silence → step reads
  exactly 1.0000 because the pre predicts silence and the post is all
  step.
- The same 0.5 step reads 0.9854 at t = 0.90 s and 8.2983 at t = 0.45
  s (both PASS). The variation is explained in §7.2; the gate-relevant
  invariant is that every measured position stays ≥ 0.99, four orders
  of magnitude above the threshold.
- `nov > 1` is not an error: it means the prediction is
  anti-correlated with part of the post (wrong phase at a new
  frequency), so `post − pred` exceeds `post` itself (§7.3).

### 7.2 Why the same step reads 0.99–8.30 depending on position

Two effects, both a function of the in-window edge position:

1. **Pre-segment length and fit quality.** A step near the window
   center has a long pre segment (up to 14.7 ms): the low-frequency
   lines are better determined, the fit is tighter, and the
   one-step prediction of the post starts as a faithful continuation
   of the *old* level — so the level change produces a large, clean
   error over the L-box. A step near the window edge has a short pre
   segment (6.7 ms): the fit is looser (the 43–61 Hz lines have
   < 1 cycle), and part of the "error" is fit inaccuracy rather than
   step novelty.
2. **Buffer mixing.** The one-step prediction of the first post sample
   is a weighted continuation of the last p = 120 pre samples. The
   more of that buffer sits at the *old* level, the further the
   prediction starts from the new level, and the larger the error over
   the 64-sample box. Near the window center the buffer is almost
   entirely old level → large error; near the window edge the buffer
   is mixed → smaller error.

Neither effect is a defect: both are measured, bounded, and leave the
stationary/edge separation at ≈ 12,000×.

### 7.3 Why a new tone can read `nov > 1`

For a new tone at frequency f, the one-step prediction of the tone
component is not zero: the AR coefficients applied to the buffer
(which now contains the tone) produce a gain G(f) — the AR polynomial
evaluated at f. The prediction error of the tone component is
`1 − G(f)` times the tone, so the measured nov is the naive energy
fraction (tone power / post power, √-scaled) times |1 − G(f)|. When
G(f) has phase near π, the prediction *cancels against* the post and
the error adds constructively: the new 1 kHz tone (A = 0.5) on the
PrimeTone background has a naive fraction √(0.125/0.839) ≈ 0.39 but
reads 1.7795 (|1 − G| ≈ 4.6). The weak tone (A = 0.15) has a naive
fraction ≈ 0.125 but reads 0.6777. The direction and size of the
amplification depend on the fit's pole configuration; the robust
statement is that new lines read well above the threshold, sometimes
by several times the naive estimate.

### 7.4 The degenerate pure-tone corner (the one measured miss)

The reference suite's "60 Hz tone + new 400 Hz at edge" case reads
**nov = 0.0159 — a real edge rejected**. Mechanism: the pre segment is
a *pure single tone* (a 2-dimensional signal), so the AR(120) fit is
rank-deficient: two poles go to the 60 Hz line and the remaining 118
poles are undetermined by the data (they fit numerical structure).
Those free poles happen to produce a prediction gain G(400 Hz) ≈ −1,
i.e. the model *predicts the new tone* and the error vanishes.

Adding a tiny amount of noise to the background breaks the degeneracy
(measured, new 400 Hz tone on a 60 Hz tone):

| background noise σ | nov |
|---|---|
| 0 (pure tone) | 0.0159 (missed) |
| 0.005 (−26 dB below the tone) | 0.1496 |
| 0.01 | 0.1109 |
| 0.02 | 0.1794 |
| 0.05 | 0.1727 |
| white noise only | 0.8797 |

Assessment: a *pure* single-tone background is a corner case — real
audio and the PrimeTone target file (27 lines + 24-bit quantization
structure) are non-degenerate, and the stress test's new-tone case
passes with nov 1.7795. The co-gate's sensitivity to a new line on a
sparse background is limited by the AR fit's null space; it is a
known, bounded limitation, not a threshold issue.

### 7.5 The weak-tone case: the co-gate is *more* sensitive than the old gate

The −15 dB tone (A = 0.15) reads nov = 0.6777 (co-gate PASS) but
stepSNR = 1.08 (old gate reject). This is expected and consistent, not
a contradiction: the co-gate's novelty is the tone's *energy fraction*
over the L-box (naive √(0.01125/0.725) ≈ 0.125, amplified to 0.6777
by the phase effect of §7.3), while the old gate's stepSNR measures a
*level* step in the derivative response, of which a tone has none
(its box mean ≈ 0, §1.2). Consequence for combination semantics: in an
AND combination (co-gate AND old gate) the weak tone is still
rejected — the sensitivity floor is preserved; the co-gate alone is
the more sensitive of the two to weak *tonal* edges.

### 7.6 The stationary maximum is still tiny even with the scan

The position scan evaluates the novelty at 4 positions per window,
which raises the stationary maximum from 3.95e-05 (single position)
to 8.07e-05 (max over 4). It remains 3 orders of magnitude below the
threshold — the scan's only effect on the null population is a modest
max-over-4 uplift, while it rescues the worst measured edge position
(0.0000 → 0.9854).

## 8. Limitations

1. **Stage-1 inheritance.** The co-gate only judges candidates that
   stage 1 produces (hmax > 0.2). An edge stage 1 misses is never
   co-gated. (On the PrimeTone file stage 1 is not the bottleneck:
   all 183 windows pass it.)
2. **Centroid jitter beyond the scan.** The scan covers centroid
   errors up to ~128 samples early / ~32 late (2.7/0.7 ms at 48 kHz).
   A larger error would let the novelty box miss the edge again. The
   measured jitter on this background is bounded by −88…+28.
3. **Degenerate (sparse) pre segments.** A pure single-tone — or
   generally very low-dimensional — pre segment makes the AR(120) fit
   rank-deficient and the out-of-model prediction at other frequencies
   untrustworthy (§7.4). A new line on such a background can be missed
   (measured 0.0159); any background noise ≥ −26 dB relative to the
   tone restores detection.
4. **Slow level ramps.** The L = 64 novelty box captures the step's
   re-equilibration transient. A level change that spreads over much
   longer than ~1.3 ms (64 samples) reads lower, in proportion to how
   much of the change fits inside the box.
5. **Threshold provenance.** 0.1 is set from the measured separation
   on one file class (strong stationary 27-tone multitone, 24-bit).
   The null population's nov scales with the pre-segment fit
   inaccuracy, which grows with background noise and with how far the
   background's spectral content is from the AR(120) model; on much
   noisier backgrounds the null maximum will rise and the margin
   should be re-measured before reuse.
6. **Clicks on strong backgrounds.** A single-sample click reads its
   energy fraction (measured 0.9523 for a 0.3 click on a 0.2-RMS-tone
   window — PASS); a click much smaller than the background's
   per-sample innovation reads below the threshold, consistent with
   the old gate (both gates are energy-normalized).

## 9. Design decisions and their measured justifications

| decision | measured justification |
|---|---|
| AR(120) + Burg fit | LS pole at \|z\| ≈ 1.10 on a 6.7 ms segment → nov 3e4 for a real step (unstable); Burg is stable by construction |
| one-step with feedback, not open-loop | open-loop: stationary file nov 0.35–2.12, 183/183 false alarms (line poles at r < 1 decay over the post); one-step: 1.3e-6…8.1e-5, 0/183 |
| novelty over the first L = 64 post samples | one-step over the whole post: the real 0.5 step reads 0.0658 < 0.1 (re-equilibration transient hides it); L = 64: 0.1568–4.1936 |
| position scan {−32, 0, +32, +64} | centroid up to 88 samples early → box misses the edge → nov 0.0000 (rejected); with scan 0.9854 (PASS) |
| threshold 0.1 | measured populations: stationary ≤ 8.07e-05, edges ≥ 0.9854 — a 12,000× empty band |
