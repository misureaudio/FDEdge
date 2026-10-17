%% replay_riesz.m - reproduce golden/results_riesz.txt via fdedge (rng(0))
% M10: thin replay driver for the riesz PoC (edge_fd_riesz_track.m).
%
% Signal generation and the fprintf blocks are the PoC's, verbatim (the
% format strings are copied line-for-line from edge_fd_riesz_track.m); the
% per-frame pipeline runs inside fdedge.tracker.Tracker (M7) and the
% calibration grids inside fdedge.calib.gridRiesz (M5).
%
% Riesz specifics: M-point PERIODIC grid x = (0:M-1)'/M (1024 points),
% symmetric edge A*|x - xe|^beta, one randn(M,1) draw per frame (RNG
% order = PoC's: rng(0), calibration is noise-free, 40 draws of
% randn(1024,1)).
%
% Pass structure (as in the PoC): the raw pass (pass 1) generates the
% signals and collects the raw estimates; the KF (pass 2) re-runs
% Tracker.update on the stored signals with the SAME xguess sequence
% (update is deterministic in (signal, xguess) -> identical raws) and the
% PoC-identical 4-state KF (kfR = [Rx 0; 0 Rb] gives the PoC's explicit
% P0 = [Rx, Rx*1e2, Rb, Rb*1e4] inside the Tracker).
%
% Frames are passed RAW to Tracker.update: the Tracker smooths internally
% (M7 contract) - pre-smoothing would double-smooth.
%
% PNG plots are skipped (visual only, not gated).
%
% Run via tests/test_replay_riesz.m (gate_replay runs this in a temp dir).

root = 'D:/Source/hermes-dir/fdedge';
addpath(fullfile(root, 'fdedge'));
addpath(fullfile(root, 'fdedge', 'fd'));
addpath(fullfile(root, 'fdedge', 'est'));
addpath(fullfile(root, 'fdedge', 'calib'));
addpath(fullfile(root, 'fdedge', 'tracker'));
addpath(fullfile(root, 'fdedge', 'signal'));

rng(0);
cfg = config('Profile','riesz');
A = cfg.A;
M = cfg.M; x = (0:M-1)'/M;   % M-point periodic grid on [0,1)
nfr = cfg.nfr;
xe    = cfg.xe0 + cfg.v*(0:nfr-1)';
betat = cfg.b0 + cfg.w*(0:nfr-1)';
alphas = cfg.alphas;
sigB = cfg.sigB;

fid = fopen('results_riesz.txt','w');
fprintf(fid, '=== Two-sided Riesz (-Delta)^{a/2} edge tracker, rng(0) ===\n');
fprintf(fid, 'x_e: 0.20 -> 0.59 (v = 0.01/frame); beta: 0.55 -> 0.85 (w = %.5f/frame)\n', cfg.w);
fprintf(fid, 'windowed-FFT fractional Laplacian on the full [0,1] domain (M = %d, periodic grid)\n', M);

% ---------- calibration: noise-free bias grids (no RNG) ----------
xes = (0.20:0.10:0.60)';
calib = gridRiesz(cfg, xes, alphas);
fprintf(fid, '\ncalibration delta_beta(xref, beta) (noise-free, in beta units):\n');
for ie = 1:numel(xes)
    line = '';
    for ib = 1:numel(alphas)
        line = [line sprintf(' %+.3e', calib.dbeta(ie,ib))];
    end
    fprintf(fid, '  xref = %.2f :%s\n', xes(ie), line);
end
fprintf(fid, '\ncalibration delta_x(xref, alpha, beta) in grid units (M = %d):\n', M);
for ia = 1:numel(alphas)
    fprintf(fid, '  alpha = %.2f :\n', alphas(ia));
    for ie = 1:numel(xes)
        line = '';
        for ib = 1:numel(alphas)
            line = [line sprintf(' %+6.1f', calib.dxg(ie,ia,ib)*M)];
        end
        fprintf(fid, '    xref = %.2f :%s\n', xes(ie), line);
    end
end

% ---------- raw per-frame estimates (pass 1) ----------
% NOTE: Tracker.update takes the RAW frame and smooths internally (M7
% contract) - do NOT pre-smooth here (double smoothing).
trk0 = Tracker(cfg, calib);
xhat_raw = zeros(nfr,1);
bhat_raw = zeros(nfr,1);
xguess = xe(1);   % initialization: known trajectory start
xguess_seq = zeros(nfr,1);
sigs = zeros(nfr, M);
for n = 1:nfr
    xguess_seq(n) = xguess;
    fe = A*abs(x - xe(n)).^betat(n);
    fn = fe + sigB*randn(M,1);
    sigs(n,:) = fn';
    [st, trk0] = trk0.update(fn, xguess);
    xhat_raw(n) = st.xraw;
    bhat_raw(n) = st.braw;
    % const-velocity guess update for the next frame (PoC: on a
    % degenerate frame the estimate is carried over and
    % xguess = xhat_raw(n) - the continue skips the update below)
    if ~isfinite(xguess)
        xguess = xhat_raw(n);
    else
        if n > 1
            xguess = xhat_raw(n) + (xhat_raw(n) - xhat_raw(n-1));
        else
            xguess = xhat_raw(n);
        end
    end
end

% ---------- Kalman (pass 2): PoC-identical 4-state KF ----------
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

fprintf(fid, '\n--- per-frame (every 5th) ---\n');
fprintf(fid, ' n   x_e_true  x_raw     x_kf      beta_true  b_raw     b_kf\n');
for n = 1:5:nfr
    fprintf(fid, '%2d  %.4f    %.4f    %.4f    %.4f     %.4f    %.4f\n', ...
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
