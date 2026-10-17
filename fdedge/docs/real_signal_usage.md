# Using fdedge on real signals — the audio case

**Status.** Written from the v1.0.0 discussion (2026-10-17). Every number
below was measured headless on base MATLAB R2024a (no toolboxes) against
the v1.0.0 library. The PoC's own signals (zero baseline, pure cusp +
white noise) are the reference case; real signals need the two
adaptations in §3.

## 1. The mapping

fdedge is a **window-local one-sided power-law cusp detector +
characterizer**: given a signal window containing an edge
`f ≈ A·(t − t_e)^β` (plus baseline and noise), it returns a sub-grid
edge position and a local exponent β. An audio transient (percussive
onset, click, attack) is exactly this model: a power-law rise riding on
a DC level + tones + noise.

| fdedge concept | audio counterpart |
|---|---|
| frame (M+1 points on `(0:M)'/M`) | analysis window of M+1 samples starting at sample index `n0` |
| one-sided edge `A·(x − x_e)^β` | transient onset `A·(t − t_e)^β` |
| exponent β | onset sharpness (β < 1 = non-differentiable cusp; smaller = sharper) |
| white noise σ | recording noise |
| Tracker across frames | sliding / overlapping-window analysis with the KF |

**Profile choice: `beta`.** Its calibration (`gridSlopePeak`) is built
on a one-sided edge `(max(x − xref, 0))^β_ref`, and its Caputo–GL
operator is one-sided — no periodicity assumption. (The Riesz path
treats the window as one period: a window cut mid-signal wraps around
and adds a P5-class boundary artifact — wrong tool for audio windows.)

## 2. Recipe (verified end-to-end)

```matlab
[audio, Fs] = audioread('clip.wav');          % base MATLAB, no toolbox
audio = audio(:,1);                           % mono
M    = min(Mmax, 8192);                       % P4 guard: M <= 8192 (186 ms @ 44.1 kHz)
n0   = T0 - round(0.4*M);                     % transient near window center (x* ~ 0.4)
w    = audio(n0:n0+M);                        % M+1 samples (grid contract, P2)
w    = w - median(w(1:round(0.25*(M+1))));    % REQUIRED for real signals: pre-edge baseline

xguess = (T0 - n0)/M;                         % coarse locator -> relative position
cfg    = config('Profile','beta','M',M);      % rebuild calibration per M (D7)
calib  = gridSlopePeak(cfg);
trk    = Tracker(cfg, calib);
[st, trk] = trk.update(w, xguess);            % frame 1 -> st.xraw / st.braw (no KF)

sample_hat = n0 + st.xraw*M;                  % back-transform
t_hat      = (sample_hat - 1)/Fs;
```

## 3. Two requirements for real signals (absent from the PoC)

### 3.1 Baseline subtraction (new finding — not in pitfalls P1–P8)

The finite-window GL does **not** kill a DC level. The window start
acts as a step edge (0 → C), and its derivative is a decaying power-law
memory tail ∝ x^(−a). Measured on the v1.0.0 library (C = 0.5, a = 0.7,
sGrid = 8, M = 1024):

| x     | 0.05  | 0.10  | 0.20  | 0.30  | 0.39  | 0.50  | 0.70  | 0.90  |
|-------|-------|-------|-------|-------|-------|-------|-------|-------|
| R     | 0.718 | 0.419 | 0.250 | 0.187 | 0.155 | 0.130 | 0.103 | 0.086 |

The tail decays as exactly x^(−0.7) (the GL of a step) and sits at
0.13–0.16 across the c1 probe window `[xguess+0.03, xguess+0.28]` —
about 20% of the edge response (0.6–0.8 at A = 0.8). It adds a slow
**negative** apparent slope to the log-log c1 fit, biasing β̂ low, and
raises the noise floor of the position metric. (The same zero-padding
mechanism (P3) also produces a large negative response at the window
end — R = −1.95 at x = 0.99 for this C — another reason the probe
windows must stay interior.)

**Fix:** subtract the pre-edge baseline before the pipeline:

```matlab
w = w - median(w(1:round(0.25*(M+1))));   % robust; uses only pre-onset samples
```

Effect (50-draw study, §4): position jitter −40%, exponent bias removed.

### 3.2 Coarse locator (xguess)

The position search is restricted to **±0.02 around `xguess`** (±20
samples at M = 1024 — the P3 boundary-restriction machinery). If the
locator misses by more than that, the estimate **pins to the search
boundary**: measured staircase — the error is constant over ±10-sample
locator ranges, jumping in ~9-sample steps. A simple energy-peak or
onset pre-detector providing ±15 samples is sufficient.

## 4. Verified performance

Single-frame study, 50 draws each, M = 1024, σ = 0.02 (the PoC's own
noise level), A = 0.8, x* = 0.390625, xguess = x* + 3 samples:

| Signal | Position err (samples) | Exponent err |
|---|---|---|
| PoC model: cusp + white noise | mean +0.01, std **2.9**, range [−11, +7] | mean −0.001, std 0.023 |
| Audio, baseline-subtracted | mean +0.8, std **3.7**, range [−7, +8] | mean −0.003, std 0.026 |
| Audio, DC left in (no baseline) | mean −0.7, std **5.2**, range [−13, +11] | mean **−0.036**, std 0.026 |

Read: **unbiased with baseline subtraction**; a few samples per frame
at the PoC's SNR. The no-HP row isolates both effects of the DC tail:
+40% position jitter and a −0.036 exponent bias (skewed distribution —
worst draw −0.088). A single bad draw can be ~5σ: one observed run gave
−18 samples and β̂ = 0.55 (clamped to the bank floor).

**The KF across frames is the intended denoiser** (that is how the PoC
made raw jitter converge). For a drifting feature, feed successive
windows to the same `Tracker`: the const-velocity KF adapts to the
drift from the raw estimates (in window-relative units a sliding window
over a static feature drifts at −hop/M). `cfg.v` is *not* used inside
the Tracker — it is the replay driver's trajectory parameter. For a
static transient only ~4–7 independent overlapping frames exist before
the edge exits the window (M = 1024, hop = 64, x* starting at 0.4), so
single-frame + a re-check window is the realistic pattern.

## 5. Constraints carried over from v1

- **M ≤ 8192** (P4 guard); rebuild the calibration for every M (D7):
  `calib = gridSlopePeak(cfg)` with `cfg.M = M`.
- **Smaller M forces smaller smoothing** (D9: `w0 > 3·sGrid/M`):
  M = 256 → sGrid ≤ 2.
- **Grid convention (P2):** M+1 samples on `x = (0:M)'/M` for the
  beta profile.
- **Window placement:** x* ∈ [0.3, 0.5] (calibration sweet spot; the
  hybrid/riesz debias grids only cover xref ∈ [0.2, 0.6]).
- **Probe window interior:** keep `xguess + w1 ≤ 1` (default
  w1 = 0.28 → xguess ≤ 0.72) — automatically satisfied when x* ≤ 0.5
  and the locator is within ±0.02.

## 6. What is missing for production audio (v2 territory)

- **No detection layer.** v1 assumes the transient is in the window
  *and* the locator is within ±0.02. Real use needs a sweep: slide a
  short detector window, rank the response energy, non-maximum
  suppress — the "multi-edge scenes" v2 item, buildable on the same
  library functions.
- **Model = cusp + white noise + flat baseline.** Tones and slow
  trends add small backgrounds (the 3 Hz tone in the study contributed
  ~0.004 to the response — negligible); strong linear trends should be
  detrended at the application level.
- **Per-frame std ~3–4 samples at the PoC's σ** (σ = 0.02, A = 0.8):
  accept it and average across windows, or improve the SNR.
- See README §"v2 candidates" for the full list.

## 7. Reproducing the studies

```matlab
M = 1024; x = (0:M)'/M; xstar = 0.390625; xguess = xstar + 3/M;
cfg = config('Profile','beta','M',M);  calib = gridSlopePeak(cfg);
Fs = 44100; n0 = 1600; Ax = 0.8; beta = 0.70; sigB = 0.02;
rng(123);                                % 456 for the no-HP row
for k = 1:50
    y = Ax*(max(x - xstar, 0)).^beta ...
      + 0.5 + 0.02*sin(2*pi*3*(n0:n0+M)'/Fs) + sigB*randn(M+1,1);
    y = y - median(y(1:round(0.25*(M+1))));   % drop this line for the no-HP row
    trk = Tracker(cfg, calib);
    [st, trk] = trk.update(y, xguess);
    % collect st.xraw, st.braw; report mean/std/range of (xraw - xstar), (braw - beta)
end
```

All studies ran headless (base MATLAB R2024a, no toolboxes) with
`rng`-seeded draws, calling the v1.0.0 flat functions
(`config`, `gridSlopePeak`, `Tracker`).
