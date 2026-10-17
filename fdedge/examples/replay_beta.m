%% replay_beta.m - reproduce golden/results_beta.txt via fdedge (rng(0))
% M9: thin replay driver for the beta PoC (edge_fd_beta_track.m).
%
% Signal generation and the fprintf blocks are the PoC's, verbatim (the
% format strings are copied line-for-line from edge_fd_beta_track.m); the
% per-frame pipeline runs inside fdedge.tracker.Tracker (M7) and the
% calibration grid inside fdedge.calib.gridSlopePeak (M5).
%
% Pass structure (as in the PoC): the raw pass (pass 1) generates the
% signals (RNG order = PoC's: rng(0), calibration is noise-free, 40 draws
% of randn(1025,1)) and collects the raw estimates; the KF (pass 2)
% re-runs Tracker.update on the stored signals with the SAME xguess
% sequence (update is deterministic in (signal, xguess) -> identical
% raws) and the PoC-identical 4-state KF (kfR = [Rx 0; 0 Rb] gives the
% PoC's explicit P0 = [Rx, Rx*1e2, Rb, Rb*1e4] inside the Tracker).
%
% Frames are passed RAW to Tracker.update: the Tracker smooths internally
% (M7 contract) - pre-smoothing would double-smooth.
%
% TABLE-SPACING NOTE (golden authority, D1): the M1 golden (the original
% PoC run, verified byte-identical to FD_EdgeDetector_PoC/results_beta.txt)
% prints the per-frame table with 2 literal spaces after %2d for the
% n = 1:5:11 rows and 3 spaces for the n = 16:5:36 rows - no single format
% string produces that, so the original script printed the table in two
% halves with different format strings. The current edge_fd_beta_track.m
% on disk is a later consolidated version ('%2d  %.4f', 2 spaces
% throughout) that does NOT reproduce the golden; the driver reproduces
% the golden's exact bytes (conditional below).
%
% PNG plots are skipped (visual only, not gated).
%
% Run via tests/test_replay_beta.m (gate_replay runs this in a temp dir).

root = 'D:/Source/hermes-dir/fdedge';
addpath(fullfile(root, 'fdedge'));
addpath(fullfile(root, 'fdedge', 'fd'));
addpath(fullfile(root, 'fdedge', 'est'));
addpath(fullfile(root, 'fdedge', 'calib'));
addpath(fullfile(root, 'fdedge', 'tracker'));
addpath(fullfile(root, 'fdedge', 'signal'));

rng(0);
cfg = config('Profile','beta');
A = cfg.A;
M = cfg.M; x = (0:M)'/M;
nfr = cfg.nfr;
xe    = cfg.xe0 + cfg.v*(0:nfr-1)';
betat = cfg.b0 + cfg.w*(0:nfr-1)';
alphas = cfg.alphas;
sigB = cfg.sigB;

fid = fopen('results_beta.txt','w');
fprintf(fid, '=== FD edge tracker with beta tracking, rng(0) ===\n');
fprintf(fid, 'x_e: 0.20 -> 0.59 (v = 0.01/frame); beta: 0.55 -> 0.85 (w = %.5f/frame)\n', cfg.w);

% ---------- calibration: slope-peak bias grid (noise-free, no RNG) ----------
dgrid = gridSlopePeak(cfg);
beta_ref = [0.55 0.65 0.75 0.85]';   % gridSlopePeak's axis (PoC)
fprintf(fid, '\ncalibration delta(alpha, beta_ref) in grid units (M = %d):\n', M);
for ia = 1:numel(alphas)
    line = '';
    for ib = 1:numel(beta_ref)
        line = [line sprintf(' %+.2e', dgrid(ia,ib)*M)];
    end
    fprintf(fid, '  alpha = %.2f :%s\n', alphas(ia), line);
end

% ---------- raw per-frame estimates (pass 1) ----------
% NOTE: Tracker.update takes the RAW frame and smooths internally (M7
% contract) - do NOT pre-smooth here (double smoothing).
trk0 = Tracker(cfg, dgrid);
xhat_raw = zeros(nfr,1);
bhat_raw = zeros(nfr,1);
xguess = xe(1);   % initialization: known trajectory start
xguess_seq = zeros(nfr,1);
sigs = zeros(nfr, M+1);
for n = 1:nfr
    xguess_seq(n) = xguess;
    fe = A*max(x - xe(n), 0).^betat(n);
    fn = fe + sigB*randn(M+1,1);
    sigs(n,:) = fn';
    [st, trk0] = trk0.update(fn, xguess);
    xhat_raw(n) = st.xraw;
    bhat_raw(n) = st.braw;
    % xguess update for the next frame (PoC: const-velocity extrapolation;
    % on a degenerate frame the estimate is carried over and
    % xguess = xhat_raw(n) - the ka/kb test below is the PoC's, re-evaluated
    % on the same (unchanged) xguess)
    ka = max(2, round((xguess + cfg.w0)*M) + 1);
    kb = min(M+1, round((xguess + cfg.w1)*M) + 1);
    if ~isfinite(xguess) || ka >= kb
        xguess = xhat_raw(n);
    else
        if isfinite(xhat_raw(n))
            if n > 1
                xguess = xhat_raw(n) + (xhat_raw(n) - xhat_raw(n-1));
            else
                xguess = xhat_raw(n);
            end
        end
    end
end

% ---------- Kalman (pass 2): PoC-identical 4-state KF ----------
Rx = var(xhat_raw - xe);
Rb = var(bhat_raw - betat);
trk = Tracker(cfg, dgrid);
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

fprintf(fid, '\n--- per-frame (every 5th) ---\n');
fprintf(fid, ' n   x_e_true  x_raw     x_kf      beta_true  b_raw     b_kf\n');
for n = 1:5:nfr
    % Golden table spacing (TABLE-SPACING NOTE): the original script
    % printed the table halves with different format strings.
    if n < 16
        nfmt = '%2d  %.4f    %.4f    %.4f    %.4f     %.4f    %.4f\n';
    else
        nfmt = '%2d   %.4f    %.4f    %.4f    %.4f     %.4f    %.4f\n';
    end
    fprintf(fid, nfmt, ...
        n, xe(n), xhat_raw(n), xhat_kf(n), betat(n), bhat_raw(n), bhat_kf(n));
end
fprintf(fid, '\n--- summary ---\n');
fprintf(fid, 'position: raw RMSE = %.3e (max %.3e) | KF RMSE = %.3e (max %.3e)\n', ...
    rmse_x_raw, max(abs(xhat_raw-xe)), rmse_x_kf, max(abs(xhat_kf-xe)));
fprintf(fid, 'exponent: raw RMSE = %.3e (max %.3e) | KF RMSE = %.3e (max %.3e)\n', ...
    rmse_b_raw, max(abs(bhat_raw-betat)), rmse_b_kf, max(abs(bhat_kf-betat)));
fprintf(fid, 'KF final velocity   v = %.5f (true 0.01000)\n', trk.kf.state(2));
fprintf(fid, 'KF final drift rate w = %.5f (true %.5f)\n', trk.kf.state(4), cfg.w);

fprintf(fid, '\nDONE\n=== end ===\n');
fclose(fid);
disp('DONE');
