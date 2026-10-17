%% demo.m - out-of-PoC generality demo (rng(42))
% M12: proves the library generalizes OUTSIDE the PoC parameters (none
% of the four PoCs ran this scenario):
%   - symmetric edge with A = 2.0 (all PoCs used A = 1.0)
%   - beta drift 0.50 -> 0.892, WIDER than the order bank
%     (0.55..0.85) - exercises the exponentBank clamp at both ends
%   - xe: 0.15 -> 0.70 (v = 0.55/49 per frame, nfr = 50) - wider than
%     the PoC's 0.20 -> 0.59; the const-velocity xguess clamp
%     [xLo, xHi] = [0.17, 0.83] engages near the start
%   - sigma = 0.02, rng(42) (the PoCs used rng(0))
%
% DEVIATION from the plan's literal config line ('v',0.015): with
% v=0.015 and nfr=50 the trajectory would end at xe = 0.885, i.e.
% within w1 = 0.28 of the domain end. There the hybrid exponent
% metric window [xref+w0, xref+w1] runs off the grid and the
% zero-padded smoothing boundary artifact (P3) corrupts the log-log
% slope - the exponent estimate collapses to the bank's lower clamp
% (measured at M12: b_raw = 0.5489 vs true 0.82-0.86 for xe > 0.75).
% That is a genuine v1 boundary limitation (documented in the README);
% the demo therefore uses the plan's stated extent xe: 0.15 -> 0.70
% (v = 0.55/49), which also matches the plan's calibration extent.
%
% The calibration is REBUILT on the new extent (D7 - the calibrator is
% parameterized, no golden involved): xes = 0.15:0.10:0.70 (6 refs),
% betas = 0.50:0.05:0.90 (9 refs) via gridHybrid.
%
% Acceptance (printed, loose on purpose - a generality smoke, not a
% golden test):
%   position KF RMSE < 3e-3
%   |bhat_KF - beta_true| < 0.05 for frames 10..50
%   no NaNs in any estimate
%
% Outputs: demo_results.txt + demo.png (via the report layer).

% Project root: hardcoded, as in the four replay drivers (this
% environment's wrapper run() reports mfilename('fullFile') as a bare
% basename, so it cannot be derived portably from here).
root = 'D:/Source/hermes-dir/fdedge';
if ~isfile(fullfile(root, 'fdedge', 'config.m'))
    error('fdedge:demo', 'project root not found at %s', root);
end
addpath(fullfile(root, 'fdedge'));
addpath(fullfile(root, 'fdedge', 'fd'));
addpath(fullfile(root, 'fdedge', 'est'));
addpath(fullfile(root, 'fdedge', 'calib'));
addpath(fullfile(root, 'fdedge', 'tracker'));
addpath(fullfile(root, 'fdedge', 'signal'));
addpath(fullfile(root, 'fdedge', 'report'));

rng(42);
cfg = config('Profile','hybrid','A',2.0,'xe0',0.15,'v',0.55/49,'b0',0.50,'w',0.008,'nfr',50);
A = cfg.A;
M = cfg.M; x = (0:M)'/M;
nfr = cfg.nfr;
xe    = cfg.xe0 + cfg.v*(0:nfr-1)';
betat = cfg.b0 + cfg.w*(0:nfr-1)';
alphas = cfg.alphas;
sigB = cfg.sigB;

% ---------- calibration on the NEW extent (D7: parameterized) ----------
xes   = (0.15:0.10:0.70)';
betas = (0.50:0.05:0.90)';
calib = gridHybrid(cfg, xes, betas);

% ---------- raw per-frame estimates (pass 1) ----------
% Tracker.update takes the RAW frame and smooths internally (M7
% contract) - do NOT pre-smooth here (double smoothing).
trk0 = Tracker(cfg, calib);
xhat_raw = zeros(nfr,1);
bhat_raw = zeros(nfr,1);
xguess = cfg.xe0;   % initialization: known trajectory start
xguess_seq = zeros(nfr,1);
sigs = zeros(nfr, M+1);
for n = 1:nfr
    xguess_seq(n) = xguess;
    fe = A*abs(x - xe(n)).^betat(n);
    fn = fe + sigB*randn(M+1,1);
    sigs(n,:) = fn';
    [st, trk0] = trk0.update(fn, xguess);
    xhat_raw(n) = st.xraw;
    bhat_raw(n) = st.braw;
    % const-velocity guess update for the next frame (hybrid: no
    % degenerate-guess branch - xguess stays finite)
    if n > 1
        xguess = xhat_raw(n) + (xhat_raw(n) - xhat_raw(n-1));
    else
        xguess = xhat_raw(n);
    end
end

% ---------- Kalman (pass 2) ----------
Rx = var(xhat_raw - xe);
Rb = var(bhat_raw - betat);
trk = Tracker(cfg, calib);
trk.kfR = [Rx 0; 0 Rb];
xhat_kf = zeros(nfr,1);
bhat_kf = zeros(nfr,1);
for n = 1:nfr
    [st, trk] = trk.update(sigs(n,:)', xguess_seq(n));
    xhat_kf(n) = st.xhat;
    bhat_kf(n) = st.bhat;
end

rmse_x_raw = sqrt(mean((xhat_raw - xe).^2));
rmse_x_kf  = sqrt(mean((xhat_kf - xe).^2));
rmse_b_raw = sqrt(mean((bhat_raw - betat).^2));
rmse_b_kf  = sqrt(mean((bhat_kf - betat).^2));
maxExpErr  = max(abs(bhat_kf(10:nfr) - betat(10:nfr)));
nNaN       = sum(~isfinite([xhat_raw; bhat_raw; xhat_kf; bhat_kf]));

% ---------- report layer: demo_results.txt + demo.png ----------
fid = fopen('demo_results.txt','w');
lines = {
'=== fdedge demo: out-of-PoC scenario (hybrid profile, rng(42)) ===';
sprintf('symmetric edge f = A| x - x_e |^beta, A = %.1f', A);
sprintf('x_e: %.2f -> %.3f (v = %.3f/frame, nfr = %d); beta: %.2f -> %.3f (w = %.3f/frame)', ...
    xe(1), xe(end), cfg.v, nfr, betat(1), betat(end), cfg.w);
sprintf('calibration rebuilt on the new extent: xes = %.2f:0.10:%.2f (%d refs), betas = %.2f:0.05:%.2f (%d refs)', ...
    xes(1), xes(end), numel(xes), betas(1), betas(end), numel(betas));
sprintf('order bank %.2f..%.2f (beta drift extends beyond it - clamp exercised)', alphas(1), alphas(end));
'';
'--- per-frame (every 5th) ---';
' n   x_e_true  x_raw     x_kf      beta_true  b_raw     b_kf';
};
for n = 1:5:nfr
    lines{end+1} = sprintf('%2d  %.4f    %.4f    %.4f    %.4f     %.4f    %.4f', ...
        n, xe(n), xhat_raw(n), xhat_kf(n), betat(n), bhat_raw(n), bhat_kf(n));
end
lines{end+1} = '';
lines{end+1} = '--- summary ---';
lines{end+1} = sprintf('position: raw RMSE = %.3e (max %.3e) | KF RMSE = %.3e (max %.3e)', ...
    rmse_x_raw, max(abs(xhat_raw-xe)), rmse_x_kf, max(abs(xhat_kf-xe)));
lines{end+1} = sprintf('exponent: raw RMSE = %.3e (max %.3e) | KF RMSE = %.3e (max %.3e)', ...
    rmse_b_raw, max(abs(bhat_raw-betat)), rmse_b_kf, max(abs(bhat_kf-betat)));
lines{end+1} = sprintf('KF final velocity   v = %.5f (true %.3f)', trk.kf.state(2), cfg.v);
lines{end+1} = sprintf('KF final drift rate w = %.5f (true %.3f)', trk.kf.state(4), cfg.w);
lines{end+1} = '';
lines{end+1} = '--- acceptance (loose smoke criteria) ---';
acc1 = rmse_x_kf < 3e-3;
acc2 = maxExpErr < 0.05;
acc3 = nNaN == 0;
lines{end+1} = sprintf('position KF RMSE = %.3e (accept < 3e-3): %s', rmse_x_kf, ...
    ternary(acc1, 'PASS', 'FAIL'));
lines{end+1} = sprintf('|bhat_KF - beta| max over frames 10..%d = %.3e (accept < 0.05): %s', ...
    nfr, maxExpErr, ternary(acc2, 'PASS', 'FAIL'));
lines{end+1} = sprintf('NaN count = %d (accept 0): %s', nNaN, ternary(acc3, 'PASS', 'FAIL'));
writeResults(fid, lines);
fclose(fid);
plotTracks(1:nfr, xe, xhat_raw, xhat_kf, betat, bhat_raw, bhat_kf, 'demo.png', ' (demo)');

fprintf('demo: position KF RMSE = %.3e | exponent max err (frames 10..%d) = %.3e | NaNs = %d\n', ...
    rmse_x_kf, nfr, maxExpErr, nNaN);
fprintf('demo acceptance: %s / %s / %s\n', ...
    ternary(acc1, 'PASS', 'FAIL'), ternary(acc2, 'PASS', 'FAIL'), ternary(acc3, 'PASS', 'FAIL'));
disp('DONE');

function s = ternary(c, a, b)
s = a;
if ~c
    s = b;
end
end
