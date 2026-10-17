# fdedge — fractional-derivative edge detector & tracker (v1)

**What it is.** `fdedge` tracks the position `x_e(n)` and the local
exponent `β(n)` of a moving power-law edge `f(x) ≈ A·|x − x_e|ᵝ` in a
noisy 1-D signal frame by frame. It does so with fractional
derivatives: a bank of Caputo (Grunwald–Letnikov) or windowed-FFT
Riesz (`(−Δ)^{a/2}`) responses is computed per frame, an edge
estimator (half-crossing, slope peak, or cusp minimum) localizes the
edge, a log-log slope or 3-parameter grid fit recovers the exponent,
noise-free calibration grids debias both, and a 2/4-state
constant-velocity Kalman filter smooths the trajectories.

**Origin.** The package is a verbatim refactor of four standalone
proof-of-concept scripts (`edge_fd_poc.m`, `edge_fd_beta_track.m`,
`edge_fd_riesz_track.m`, `edge_fd_hybrid_track.m`, in
`D:\Source\hermes-dir\FD_EdgeDetector_PoC\`, kept read-only). Every
moved expression was copied line-for-line (rename + package prefix
only); the four `examples/replay_*.m` drivers reproduce the PoCs'
published result files through the library and are gated against byte
snapshots in `golden/` (see *Verification* below).

**Requirements.** Base MATLAB R2024a, no toolboxes (`matlab.unittest`
is base MATLAB). Deterministic: replays start with `rng(0)` and consume
the RNG in the exact call order of the PoCs (one `randn` per frame,
same sizes) — replay-critical.

**Real-signal usage** (e.g. audio windows):
[docs/real_signal_usage.md](docs/real_signal_usage.md) — mapping,
verified recipe, the DC-baseline requirement, 50-draw performance.
**Sliding-window sweep** (finding edges in a long wav/flac file):
[docs/sliding_window_sweep.md](docs/sliding_window_sweep.md) —
two-stage design (step-height scan + per-candidate v1 pipeline),
verified 8/8 on 60 s of audio, ~35× realtime; reference
implementation [examples/sweep_demo.m](examples/sweep_demo.m).

---

## Quickstart

The shape of `examples/demo.m` (out-of-PoC scenario, hybrid profile):

```matlab
addpath('D:/Source/hermes-dir/fdedge/fdedge');  % + its fd/, est/, calib/,
addpath('D:/Source/hermes-dir/fdedge/fdedge/fd');   % tracker/, signal/, report/
addpath('D:/Source/hermes-dir/fdedge/fdedge/est');
addpath('D:/Source/hermes-dir/fdedge/fdedge/calib');
addpath('D:/Source/hermes-dir/fdedge/fdedge/tracker');
addpath('D:/Source/hermes-dir/fdedge/fdedge/signal');

rng(42);
cfg   = config('Profile','hybrid','A',2.0,'xe0',0.15,'v',0.55/49,'b0',0.50,'w',0.008,'nfr',50);
calib = gridHybrid(cfg, (0.15:0.10:0.70)', (0.50:0.05:0.90)');  % rebuild for your extent
trk   = Tracker(cfg, calib);
xguess = cfg.xe0;                                   % known trajectory start
for n = 1:cfg.nfr
    frame = cfg.A*abs(x - xe(n)).^betat(n) + cfg.sigB*randn(cfg.M+1,1);  % RAW 1025-pt signal
    [st, trk] = trk.update(frame, xguess);          % st: xhat, bhat, vhat, ohat (+xraw, braw)
    xguess = st.xhat + (st.xhat - xhat_prev);       % const-velocity prediction for next frame
end
```

Per-profile signal contract (enforced, error `fdedge:gridLen`):
`pocB`/`beta`/`hybrid` take 1025 points on `x = (0:1024)'/1024`;
`riesz` takes 1024 points on `xr = (0:1023)'/1024` (periodic). The
Tracker smooths internally — pass raw frames.

## API (v1 surface)

Flat function names in plain directories (see pitfall P7 — headless
MATLAB in this environment never resolves `+package` names).

| Function | Purpose |
|---|---|
| `config('Profile', p, 'field', value, ...)` | config factory; PoC-matching defaults per profile + `winScale` validation (D9) |
| `caputoGL(t, f, a)` / `caputoL1(t, f, a)` | Caputo derivative: O(N log N) GL convolution / product integration (α∈(0,1), P1 guard) |
| `rieszLaplacian(f, a \| alphas)` | windowed-FFT `(−Δ)^{a/2}`; bank → M×numel(alphas); M≤8192 (P4 guard) |
| `smoothGauss(f, sGrid)` | Gaussian pre-smoothing (`conv(...,'same')`; P3 boundary artifact) |
| `derivativeBank(t, f, alphas, 'Method', m)` | facade: `'gl'` \| `'l1'` \| `'riesz'` |
| `cuspMin(R, xg, xlo, xhi)` / `slopePeak(S, xg, xlo, xhi)` / `halfCross(R, xg, level)` | position estimators (parabolic sub-grid; window-restricted) |
| `loglogSlope(y, d)` | `polyfit(log d, log max(y,1e-6), 1)` slope |
| `memorySubtract(fs, xref, M)` | left-side truncation for the hybrid exponent (P6) |
| `exponentBank(c1, alphas)` | zero of the c1-vs-α line + bank clamp → `[bhat, iaStar]` |
| `rieszGfit(y, d)` | 3-parameter (C, g, C0) grid search → `[ghat, ok]` (P5) |
| `flatnessRatio(R, alphas, N, w, iBeta)` | Part A one-shot: endpoint ratio, plateau mean, best order |
| `interp2d(G, rAx, cAx, qr, qc)` | bilinear lookup with clipping |
| `scalarHalfCross(cfg)` / `gridSlopePeak(cfg)` / `gridRiesz(cfg, xes, betas)` / `gridHybrid(cfg, xes, betas)` | calibration builders (pocB / beta / riesz / hybrid) |
| `KalmanCV(n, qV, qW, Rmat)` | 2/4-state constant-velocity KF; `init(x0, P0)`, `predict()`, `update(z)` |
| `Tracker(cfg, calib)` | per-frame orchestrator; `update(signal, xguess)` → `xhat, bhat, vhat, ohat` (+ raw) |
| `analyzeMatchedOrder(cfg, t, f)` | Part A one-shot (L1 bank → flatness ratio) |
| `writeResults(fid, lines)` / `plotTracks(...)` | report layer (demo output; replays keep verbatim PoC fprintf blocks) |

## Verification

Run the whole suite headless: `tests/run_tests.m` → final line
`TESTS_PASSED` (77 tests: operator references, estimator
sub-grid/robustness, calibration invariances, KF invariants, per-profile
Tracker frames, and the four golden replays below).

**The two-tier gate.** Each replay driver writes its results file; the
gate (a) *stage 1 — PRIMARY, portable*: parses every printed number
from golden and replay and requires agreement within **2× the number's
last printed digit** (integers exact) — tolerates cross-platform
last-ulp flips, catches any real numerical regression; (b) *stage 2 —
SECONDARY, this environment*: raw byte `isequal`. If stage 1 passes
but stage 2 fails, the gate triages: strip CR from both and re-compare
(pure line-ending difference → environmental); any remaining byte diff
is logged in `tests/replay_exceptions.md` for review. A genuine
numeric difference always fails stage 1 and is red regardless.

| Replay | Golden | len (B) | Σ uint8 | Stage 1 | Stage 2 |
|---|---|---|---|---|---|
| `replay_poc.m` | `results_poc.txt` | 1475 | 104092 | PASS | CRLF triage (golden CRLF, replay LF; LF-normalized bytes identical) |
| `replay_beta.m` | `results_beta.txt` | 1479 | 83365 | PASS | byte-identical |
| `replay_riesz.m` | `results_riesz.txt` | 4167 | 210800 | PASS | byte-identical |
| `replay_hybrid.m` | `results_hybrid.txt` | 4271 | 219752 | PASS | byte-identical |

`tests/golden_manifest.txt` holds the len/Σ pairs (machine-readable);
`test_manifest` re-runs all four replays and checks both the manifest
and the two-tier gate. `tests/replay_exceptions.md` is absent: the only
stage-2 deviation (poc) is a pure line-ending difference, which the
gate classifies as environmental without logging.

## Pitfalls (each has a test)

- **P1** The L1 scheme degenerates at α=1 (silent zeros) → `caputoL1` rejects α≥1.
- **P2** Grids are normalized `(0:M)'/M`; a missing `/M` silently moves the edge to the boundary.
- **P3** Zero-padded Gaussian smoothing creates spurious boundary edges → all extremum searches are window-restricted (tests place boundary spikes and verify rejection).
- **P4** k-vector orientation trap: column `fft(f) .*` row `abs(k).^a` → M×M matrix → OOM; k is a column, entry points coerce `f=f(:)`.
- **P5** Windowed Riesz of a non-decaying signal = cusp + curved background; the 3-param gfit runs away to the grid edge under noise near boundaries (keep the ±0.45 range, one-sided probes, no symmetric average, no quadratic background).
- **P6** Memory subtraction `D^a f − D^a f_left = D^a g` removes the left-side memory exactly, with no added noise; robustness envelope tested at ±2–3 grid-point position error.
- **P7** (this environment) Headless/batch MATLAB never resolves `+package` names → plain directories + flat names; any future `+package` design must be re-tested headless first.
- **P8** The windowed-FFT fractional operator's `|k|^a` multiplier has a cusp at k=0, so the period-L operator converges to the continuum operator only algebraically, ~L^{−(1+a)}; refining dx at fixed L does nothing.

## Calibration sensitivity (D7)

The calibration grids are noise-free bias maps of the *pipeline*, not
of the signal: they are **A-invariant** (tested to ≤1e-12 — the
estimators are scale-invariant in A) but **sensitive to σ_smooth
(`sGrid`), M, and the probe windows (`w0`/`w1`)** — those change the
operator, hence the bias. Rebuild the grid whenever any of them
changes: the calibrators take the parameters (`gridHybrid(cfg, xes,
betas)` etc.) and recompute on the fly — `examples/demo.m` rebuilds
the grid on a new extent (6×9 refs) as its first act.

## Known limitations

- **Riesz gfit runaway near boundaries under noise (P5):** with the
  edge close to the domain end, the background fit can run away to the
  grid edge; the flatness test reports `ok=false`, but the estimate is
  then untrustworthy.
- **Exponent KF lag on drift:** the KF is a constant-velocity
  smoother; under sustained β drift it lags the raw (unbiased)
  estimates by roughly drift-rate × lag.
- **M ≤ 8192:** the sandbox cannot FFT beyond ~8192 points without OOM;
  `rieszLaplacian` refuses M>8192 explicitly. Defaults stay at M=1024.
- **Zero-padded smoothing boundary artifact (P3):** extremum searches
  are window-restricted to mitigate it, but the *hybrid exponent
  metric window* `[xref+w0, xref+w1]` must also stay inside the
  domain: with xref > 1−w1 (= 0.72 at the defaults) the window runs
  off the grid and the boundary artifact corrupts the log-log slope —
  the exponent estimate collapses to the bank's lower clamp
  (measured during M12 demo development: xe > 0.75 → b_raw ≈ 0.549
  vs true 0.82–0.86). Keep the tracked trajectory within
  `[xLo, 1−w1]` until v2 boundary handling lands.
- **Algebraic windowed-FFT convergence (P8):** Riesz results carry a
  ~L^{−(1+a)} discretization error that no dx refinement removes at
  fixed L.
- **DC level in real signals:** the finite-window GL does not remove a
  DC offset — the window-start step produces a decaying memory tail
  (∝ x^{−a}) that inflates position jitter by ~40% and biases β̂ low by
  ~0.036 (50-draw study). Baseline-subtract the window (pre-edge
  median) before the pipeline — see docs/real_signal_usage.md.

## v2 candidates (explicitly out of v1 scope)

- `'replicate'` boundary padding (D12) — would remove the P3-class
  boundary artifacts above, but changes the smoothed signal and
  therefore **invalidates all four goldens**: a full recalibration and
  new golden snapshots would be required (the cost is stated here so
  the decision is informed).
- Multi-edge scenes (per-edge tracking, edge splitting/merging).
- 2-D edges.
- Adaptive order banks (per-frame/per-region α selection).
- `fdedge-py`: Python port (D8).

## Layout

```
fdedge/
├── fdedge/            % the library: config.m, analyzeMatchedOrder.m,
│   ├── fd/            %   fd/ (operators), est/ (estimators), calib/
│   ├── est/           %   (calibration), tracker/ (Tracker, KalmanCV),
│   ├── calib/         %   signal/ (synthEdge), report/ (writeResults,
│   ├── tracker/       %   plotTracks)
│   ├── signal/
│   └── report/
├── tests/             % run_tests.m, gate_replay.m, gate_files.m,
│                      % golden_manifest.txt, test_*.m (77 tests)
├── golden/            % byte-snapshots of the 4 PoC result files
├── examples/          % replay_{poc,beta,riesz,hybrid}.m, demo.m, sweep_demo.m
├── docs/              % real_signal_usage.md, sliding_window_sweep.md
└── README.md
```restricted) |
| `loglogSlope(y, d)` | `polyfit(log d, log max(y,1e-6), 1)` slope |
| `memorySubtract(fs, xref, M)` | left-side truncation for the hybrid exponent (P6) |
| `exponentBank(c1, alphas)` | zero of the c1-vs-α line + bank clamp → `[bhat, iaStar]` |
| `rieszGfit(y, d)` | 3-parameter (C, g, C0) grid search → `[ghat, ok]` (P5) |
| `flatnessRatio(R, alphas, N, w, iBeta)` | Part A one-shot: endpoint ratio, plateau mean, best order |
| `interp2d(G, rAx, cAx, qr, qc)` | bilinear lookup with clipping |
| `scalarHalfCross(cfg)` / `gridSlopePeak(cfg)` / `gridRiesz(cfg, xes, betas)` / `gridHybrid(cfg, xes, betas)` | calibration builders (pocB / beta / riesz / hybrid) |
| `KalmanCV(n, qV, qW, Rmat)` | 2/4-state constant-velocity KF; `init(x0, P0)`, `predict()`, `update(z)` |
| `Tracker(cfg, calib)` | per-frame orchestrator; `update(signal, xguess)` → `xhat, bhat, vhat, ohat` (+ raw) |
| `analyzeMatchedOrder(cfg, t, f)` | Part A one-shot (L1 bank → flatness ratio) |
| `writeResults(fid, lines)` / `plotTracks(...)` | report layer (demo output; replays keep verbatim PoC fprintf blocks) |

## Verification

Run the whole suite headless: `tests/run_tests.m` → final line
`TESTS_PASSED` (77 tests: operator references, estimator
sub-grid/robustness, calibration invariances, KF invariants, per-profile
Tracker frames, and the four golden replays below).

**The two-tier gate.** Each replay driver writes its results file; the
gate (a) *stage 1 — PRIMARY, portable*: parses every printed number
from golden and replay and requires agreement within **2× the number's
last printed digit** (integers exact) — tolerates cross-platform
last-ulp flips, catches any real numerical regression; (b) *stage 2 —
SECONDARY, this environment*: raw byte `isequal`. If stage 1 passes
but stage 2 fails, the gate triages: strip CR from both and re-compare
(pure line-ending difference → environmental); any remaining byte diff
is logged in `tests/replay_exceptions.md` for review. A genuine
numeric difference always fails stage 1 and is red regardless.

| Replay | Golden | len (B) | Σ uint8 | Stage 1 | Stage 2 |
|---|---|---|---|---|---|
| `replay_poc.m` | `results_poc.txt` | 1475 | 104092 | PASS | CRLF triage (golden CRLF, replay LF; LF-normalized bytes identical) |
| `replay_beta.m` | `results_beta.txt` | 1479 | 83365 | PASS | byte-identical |
| `replay_riesz.m` | `results_riesz.txt` | 4167 | 210800 | PASS | byte-identical |
| `replay_hybrid.m` | `results_hybrid.txt` | 4271 | 219752 | PASS | byte-identical |

`tests/golden_manifest.txt` holds the len/Σ pairs (machine-readable);
`test_manifest` re-runs all four replays and checks both the manifest
and the two-tier gate. `tests/replay_exceptions.md` is absent: the only
stage-2 deviation (poc) is a pure line-ending difference, which the
gate classifies as environmental without logging.

## Pitfalls (each has a test)

- **P1** The L1 scheme degenerates at α=1 (silent zeros) → `caputoL1` rejects α≥1.
- **P2** Grids are normalized `(0:M)'/M`; a missing `/M` silently moves the edge to the boundary.
- **P3** Zero-padded Gaussian smoothing creates spurious boundary edges → all extremum searches are window-restricted (tests place boundary spikes and verify rejection).
- **P4** k-vector orientation trap: column `fft(f) .*` row `abs(k).^a` → M×M matrix → OOM; k is a column, entry points coerce `f=f(:)`.
- **P5** Windowed Riesz of a non-decaying signal = cusp + curved background; the 3-param gfit runs away to the grid edge under noise near boundaries (keep the ±0.45 range, one-sided probes, no symmetric average, no quadratic background).
- **P6** Memory subtraction `D^a f − D^a f_left = D^a g` removes the left-side memory exactly, with no added noise; robustness envelope tested at ±2–3 grid-point position error.
- **P7** (this environment) Headless/batch MATLAB never resolves `+package` names → plain directories + flat names; any future `+package` design must be re-tested headless first.
- **P8** The windowed-FFT fractional operator's `|k|^a` multiplier has a cusp at k=0, so the period-L operator converges to the continuum operator only algebraically, ~L^{−(1+a)}; refining dx at fixed L does nothing.

## Calibration sensitivity (D7)

The calibration grids are noise-free bias maps of the *pipeline*, not
of the signal: they are **A-invariant** (tested to ≤1e-12 — the
estimators are scale-invariant in A) but **sensitive to σ_smooth
(`sGrid`), M, and the probe windows (`w0`/`w1`)** — those change the
operator, hence the bias. Rebuild the grid whenever any of them
changes: the calibrators take the parameters (`gridHybrid(cfg, xes,
betas)` etc.) and recompute on the fly — `examples/demo.m` rebuilds
the grid on a new extent (6×9 refs) as its first act.

## Known limitations

- **Riesz gfit runaway near boundaries under noise (P5):** with the
  edge close to the domain end, the background fit can run away to the
  grid edge; the flatness test reports `ok=false`, but the estimate is
  then untrustworthy.
- **Exponent KF lag on drift:** the KF is a constant-velocity
  smoother; under sustained β drift it lags the raw (unbiased)
  estimates by roughly drift-rate × lag.
- **M ≤ 8192:** the sandbox cannot FFT beyond ~8192 points without OOM;
  `rieszLaplacian` refuses M>8192 explicitly. Defaults stay at M=1024.
- **Zero-padded smoothing boundary artifact (P3):** extremum searches
  are window-restricted to mitigate it, but the *hybrid exponent
  metric window* `[xref+w0, xref+w1]` must also stay inside the
  domain: with xref > 1−w1 (= 0.72 at the defaults) the window runs
  off the grid and the boundary artifact corrupts the log-log slope —
  the exponent estimate collapses to the bank's lower clamp
  (measured during M12 demo development: xe > 0.75 → b_raw ≈ 0.549
  vs true 0.82–0.86). Keep the tracked trajectory within
  `[xLo, 1−w1]` until v2 boundary handling lands.
- **Algebraic windowed-FFT convergence (P8):** Riesz results carry a
  ~L^{−(1+a)} discretization error that no dx refinement removes at
  fixed L.
- **DC level in real signals:** the finite-window GL does not remove a
  DC offset — the window-start step produces a decaying memory tail
  (∝ x^{−a}) that inflates position jitter by ~40% and biases β̂ low by
  ~0.036 (50-draw study). Baseline-subtract the window (pre-edge
  median) before the pipeline — see docs/real_signal_usage.md.

## v2 candidates (explicitly out of v1 scope)

- `'replicate'` boundary padding (D12) — would remove the P3-class
  boundary artifacts above, but changes the smoothed signal and
  therefore **invalidates all four goldens**: a full recalibration and
  new golden snapshots would be required (the cost is stated here so
  the decision is informed).
- Multi-edge scenes (per-edge tracking, edge splitting/merging).
- 2-D edges.
- Adaptive order banks (per-frame/per-region α selection).
- `fdedge-py`: Python port (D8).

## Layout

```
fdedge/
├── fdedge/            % the library: config.m, analyzeMatchedOrder.m,
│   ├── fd/            %   fd/ (operators), est/ (estimators), calib/
│   ├── est/           %   (calibration), tracker/ (Tracker, KalmanCV),
│   ├── calib/         %   signal/ (synthEdge), report/ (writeResults,
│   ├── tracker/       %   plotTracks)
│   ├── signal/
│   └── report/
├── tests/             % run_tests.m, gate_replay.m, gate_files.m,
│                      % golden_manifest.txt, test_*.m (77 tests)
├── golden/            % byte-snapshots of the 4 PoC result files
├── examples/          % replay_{poc,beta,riesz,hybrid}.m, demo.m, sweep_demo.m
├── docs/              % real_signal_usage.md, sliding_window_sweep.md
└── README.md
```
