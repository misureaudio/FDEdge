%% replay_poc.m - reproduce golden/results_poc.txt via fdedge (rng(0))
% M8: thin replay driver for the original PoC (edge_fd_poc.m).
%
% Signal generation and the fprintf blocks are the PoC's, verbatim (the
% format strings are copied line-for-line from edge_fd_poc.m); the
% per-frame pipeline runs inside fdedge.tracker.Tracker (M7) and the
% Part A one-shot inside fdedge.analyzeMatchedOrder (M7).
%
% Part A keeps the PoC's NON-normalized grid t = (0:N)' with A*t.^beta
% generated directly (D4: the PoC's Part A has no shifted edge) and is
% pre-smoothed here (the PoC smooths before the bank; analyzeMatchedOrder
% does not smooth).
%
% Part B KF: the PoC builds Rm = var(raw_db - xe) offline from the true
% trajectory (which it knows); this driver does a raw pre-pass (no RNG in
% Tracker.update) to collect the debiased raws, sets trk.kfR = [Rm], then
% replays the stored signals through the PoC-identical 2-state KF.
% Part B frames are passed RAW to Tracker.update: the Tracker smooths
% internally (M7 contract) - pre-smoothing would double-smooth.
%
% PNG plots are skipped (visual only, not gated).
%
% Run via tests/test_replay_poc.m (gate_replay runs this in a temp dir).

root = 'D:/Source/hermes-dir/fdedge';
addpath(fullfile(root, 'fdedge'));
addpath(fullfile(root, 'fdedge', 'fd'));
addpath(fullfile(root, 'fdedge', 'est'));
addpath(fullfile(root, 'fdedge', 'calib'));
addpath(fullfile(root, 'fdedge', 'tracker'));
addpath(fullfile(root, 'fdedge', 'signal'));

rng(0);
cfg = config('Profile','poc');
A = cfg.A; beta = cfg.beta;
fid = fopen('results_poc.txt','w');
fprintf(fid, '=== FD edge-detector PoC, rng(0), A=%.1f, beta=%.2f ===\n', A, beta);

%% ================= PART A : plateau at the matched order =============
N = cfg.N; t = (0:N)';
g  = A*(t.^beta);
gn = g + cfg.sigA*randn(N+1,1);
gs = smoothGauss(gn, cfg.sGrid);
res = analyzeMatchedOrder(cfg, t, gs);
exact_plat = res.exactPlat;
plat_mean  = res.platMean;
alphas = cfg.alphas;
fprintf(fid, '\n--- Part A (noise sigma = %.3f, smoothing %d grid pts) ---\n', cfg.sigA, cfg.sGrid);
fprintf(fid, 'exact plateau  D^beta[A t^beta] = A*Gamma(1.7) = %.6f\n', exact_plat);
fprintf(fid, 'measured plateau mean (alpha = beta) = %.6f\n', plat_mean);
fprintf(fid, 'endpoint ratio R(0.95)/R(0.05)   log-log slope (est. beta-alpha):\n');
for j = 1:numel(alphas)
    fprintf(fid, '  alpha = %.2f : ratio = %.4f   loglog-slope = %+.3f\n', alphas(j), res.ratio(j), res.c1(j));
end
fprintf(fid, 'best-matching order (min |ratio-1|): alpha = %.2f  (true beta = %.2f)\n', alphas(res.iBest), beta);
[~, iB] = min(abs(alphas - beta));
fprintf(fid, 'edge-exponent estimate at alpha=beta: beta_hat = beta - c1 = %.3f\n', beta - res.c1(iB));

%% ================= PART B : tracking a moving edge ====================
M = cfg.M; x = (0:M)'/M;
nfr = cfg.nfr;
xe = 0.25 + 0.01*(0:nfr-1)';

calib = scalarHalfCross(cfg);
% Raw pre-pass: generate the noisy frames (RNG order = PoC's) and collect
% the debiased raw estimates (st.xraw = raw - delta, as in the PoC).
% NOTE: Tracker.update takes the RAW frame and smooths internally (M7
% contract) - do NOT pre-smooth here (double smoothing).
trk0 = Tracker(cfg, calib);
raw_db = zeros(nfr, 1);
sigs = zeros(nfr, M+1);
for n = 1:nfr
    fe = A*max(x - xe(n), 0).^beta;
    fn = fe + cfg.sigB*randn(M+1,1);
    sigs(n,:) = fn';
    [st, trk0] = trk0.update(fn, NaN);
    raw_db(n) = st.xraw;
end
rmse_raw = sqrt(mean((raw_db - xe).^2));
Rm = var(raw_db - xe);
% Main pass: PoC-identical 2-state KF (Rm from the raw variances).
trk = Tracker(cfg, calib);
trk.kfR = [Rm];
kf = zeros(nfr, 1);
for n = 1:nfr
    [st, trk] = trk.update(sigs(n,:)', NaN);
    kf(n) = st.xhat;
end
rmse_kf = sqrt(mean((kf - xe).^2));

fprintf(fid, '\n--- Part B (noise sigma = %.3f, v = 0.01/frame, %d frames) ---\n', cfg.sigB, nfr);
fprintf(fid, 'estimator bias (noise-free calibration)  delta = %.2e\n', calib);
fprintf(fid, 'raw (debiased)  RMSE = %.3e   (max |err| = %.3e)\n', rmse_raw, max(abs(raw_db-xe)));
fprintf(fid, 'Kalman-filtered RMSE = %.3e   (max |err| = %.3e)\n', rmse_kf, max(abs(kf-xe)));
fprintf(fid, 'Kalman improvement factor = %.2fx\n', rmse_raw/rmse_kf);
fprintf(fid, 'Kalman final velocity estimate = %.5f (true 0.01000)\n', trk.kf.state(2));

fprintf(fid, '\nDONE\n=== end ===\n');
fclose(fid);
disp('DONE');
