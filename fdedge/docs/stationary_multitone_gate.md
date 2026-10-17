# Stationary-multitone stress test: the stepSNR gate false-positive mode

**Status.** Written 2026-10-17. Companion to
`docs/sliding_window_sweep.md` (the sweep design) and
`docs/real_signal_usage.md` (single known-edge window). All numbers
measured headless on base MATLAB R2024a against the v1.0.0 library,
plus the Python spectral characterization in
`D:\Source\hermes-dir\.venv`.

Test file: `PrimeToneNewPhi[Nch 1][Frac 03][LFL 043][00020-21600]
[Nfr 027]-[48000-24]-[1s]-[-0.25dB]-[049062]-[ 8.5905]-[RD].flac`
— a **strong stationary 27-tone multitone** (no transients anywhere
in the 1 s record). It is the worst case for a transient-edge
detector: every window contains full-strength signal.

## 1. File characterization (measured)

- FLAC PCM_24, 48 kHz, mono, 48000 frames (1.0 s), 117512 bytes.
- Peak 0.971628 (−0.250 dBFS, matches the filename tag),
  RMS −8.840 dBFS, **crest factor 8.5905 dB** (matches the filename
  tag exactly). Global peak at t = 298.188 ms.
- **27 tones** (matches the `Nfr 027` tag) at **prime-Hz
  frequencies**: 43, 53, 61, 83, 107, 127, 167, 211, 251, 337, 421,
  509, 683, 853, 1021, 1361, 1699, 2039, 2731, 3413, 4093, 5449,
  6827, 8191, 10909, 13649, 16381. A 1 s record resolves frequency to
  ±0.5–1 Hz (the low tones have only 43–127 cycles in the record);
  all measured peaks sit within that band of the stated primes
  (e.g. 43.018, 61.112, 83.004, 167.015, 16381.004). Sub-Hz
  "precision" below ~2 kHz is meaningless on this record.
- Amplitudes follow the **pink law** A ∝ f^(−0.5) (fit exponent
  −0.4985; ≤ ~1 dB scatter for the 26 integer-verified tones):
  constant power per 1/3 octave — the `Frac 03` + "pink noise
  frequency response" tags.
- Spacing is a **4/3, 5/4, 6/5 ratio cycle** (product 2 → one octave
  per 3 tones = 1/3-octave density) snapped to primes: exact from
  251 Hz upward (deviation < 0.5%), perturbed below (prime sparsity).
  The top tone 16381 stops because the next 4/3 step (≈ 21950)
  exceeds the 21600 upper bound; `LFL 043` is the bottom tone.
- **Phases (Phi):** 27 distinct phases (deg, in tone order):
  16.2, −158.8, 8.5, 99.9, 156.3, 147.4, 76.5, 63.2, 43.7, 133.1,
  13.9, 77.2, 162.6, 143.2, 150.7, 107.2, 175.7, 51.9, 157.4, 126.6,
  167.3, 98.7, 39.8, 95.4, 155.4, 65.4, 18.8. They are optimized to
  minimize crest factor: measured crest 8.5905 dB vs 500 random phase
  draws (min 9.0466, p1 9.3867, median 10.7574, max 12.3386) —
  ≥ 0.46 dB below the best random draw.
- Dominant envelope beats (measured from |y|): 10, 18, 8, 30, 74, 64,
  44, 22, 40, 54 Hz = differences of the low primes (53−43, 61−43,
  61−53, 83−53, 127−53, 107−43, 127−83, 83−61, 83−43, 107−53).

## 2. Sweep result on the file (M = 1024, H = 256, α = 0.70)

183 windows; stage 1 = 0.193 s, stage 2 = 0.63 s (per 1 s of audio).

| stage | result |
|---|---|
| 1 (hmax > 0.2) | **all 183 windows pass** — hmax median 4.352, q99 7.553, max 7.580. The step-height floor of a strong stationary tone is ~20–38× the noise-floor threshold (0.2). Stage 1 has no selectivity here by design (it is a noise-floor gate). |
| 2 (stepSNR > 3, posd ≤ 0.0201) | **24 windows pass stepSNR > 3 (20 with the position co-gate)** — stepSNR median over all candidates 1.045, **max 10.496**. Top 5: 10.496, 9.188, 7.481, 6.714, 6.299 — **all with braw = 0.550 (the bank floor)**. |

So the v1 gate **false-alarms on 13% of windows of a strong
stationary multitone** — a 10-minute file containing 1 minute of
drone/test-tone content would yield ≈ 1400 false edges. This is the
gate's verified failure mode on stationary tonal signal; the verified
60 s onset run (noise background, 0.02 σ) had max noise stepSNR 1.4
and 0 false alarms, so the two regimes are distinct.

## 3. Mechanism (measured on the top false-pass windows)

The stepSNR gate's premise is "the pre-edge reference
`std(R2(1:64))` is the noise floor." For a **stationary low tone**
(43–127 Hz = 1–2.6 cycles per 21.3 ms window) the GL response is a
slow oscillation, and the first-64-sample box can sit on its nearly
flat part while the post-edge box catches full amplitude:

| window (t) | raw pre-box std | raw post-box std | Δmean | resp pre std | resp post mean | stepSNR |
|---|---|---|---|---|---|---|
| 123 (0.651 s) | 0.0567 | 0.3018 | −0.6612 | 0.1170 | 1.2278 | 10.50 |
| 116 (0.613 s) | 0.0569 | 0.2752 | +0.7944 | 0.3188 | 2.9293 | 9.19 |
| 50 (0.261 s) | 0.0686 | 0.2968 | −0.0178 | 0.3746 | 2.8027 | 7.48 |

These are **envelope crests** of the phase-optimized multitone
(quiet phase → constructive beat peak): at the 21 ms window scale a
quiet→crest transition is a *real local level change*.

## 4. Candidate fixes tested (measured)

**(a) Raw level-change gate** L = |Δmean| / std: the 24 false passes
have L max 2.545, median 1.047 — and with pre-box normalization the
crest windows reach 11.7σ, comparable to a real step (25σ). A quiet→
crest is locally indistinguishable from a step in the level domain.
**Does not fix it.**

**(b) Box correlation** rbox = corr(pre-box, post-box) (64 samples
each): references — real step + noise **0.052**, clean step NaN
(silent pre-box), 60 Hz tone **0.947**, 1 kHz tone **0.416**, 60 Hz +
click at edge **0.947**. On the file: the 24 false passes have rbox
max **0.782**, median −0.002, min −0.894; `rbox < 0.5` rejects 20 of
24, but 4 crest windows (0.5–0.78) overlap the stationary 1 kHz tone
(0.416) — **no single rbox threshold separates them**. Partial
filter, not a fix.

**(c) Off-line spectral power ratio** (pre/post half split, line
mask from the full-window FFT): **rejected** — a step's energy is DC +
the boundary discontinuity; per-half FFTs of the smooth segments
carry none of it (measured: step reads 0.93, tones 1e17 with an
unmasked DC bin). The half-split is blind to step edges by
construction.

## 5. What actually separates the cases

Across every tested case, the separating property is **spectral
novelty**: *does the post-edge region contain spectral content the
pre-edge region does not?*

| case | new content at edge |
|---|---|
| step / level edge | DC + broadband 1/f tail — novel |
| new tone onset | a line present post, absent pre — novel |
| click / burst | broadband burst — novel |
| stationary multitone crest | **none** — the same 27 lines, re-phased |

The crest is the *only* tested case with no spectral novelty, and it
is the *only* false-positive source measured so far. A complete fix
is therefore a **per-window line-set co-gate**: full-window FFT line
set (local maxima > −30 dB, ±2-bin masks), then per-line and
off-line post/pre power ratios; reject the candidate if no line (and
no off-line band) shows a ratio above a small threshold (≈ 3–5).
Not yet implemented; the mechanism and the negative results above
are the design basis.

## 6. File-start edge (padded run)

Prepending 512 zeros (edge at x = 0.5, single-window pipeline):
xraw = 0.5276 → **t = +0.589 ms** vs true 0 (refine bias for a fast
edge, consistent with the M13 sweep), braw = 0.8364, Ahat = 1.6689,
stepSNR = Inf (pre-edge std exactly 0 — pure silence; the gate
expression needs an eps guard in `std` before use in production).
