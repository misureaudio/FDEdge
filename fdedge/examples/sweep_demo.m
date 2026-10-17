%% SWEEP_DEMO — sliding-window sweep for fdedge transient-edge detection
%
%  Drives the v1 fdedge pipeline over a long signal (e.g. an audio file)
%  to find, localize and characterize transient-related edges (onsets).
%
%  Two stages (both verified — see docs/sliding_window_sweep.md):
%
%    Stage 1  coarse scan: one Caputo-GL evaluation per window.
%             Detection statistic = STEP-HEIGHT profile (boxcar difference
%             of the GL response across a trial edge), localized by the
%             centroid of the positive profile. Threshold + NMS + dedup.
%    Stage 2  per candidate: full v1 pipeline (Tracker = alpha-bank +
%             sub-sample refine) + two-tier gate (stepSNR > 3, position
%             sanity <= refine search width).
%
%  Usage with a real file:
%      [y, Fs] = audioread('clip.wav');
%      y = y(:,1);                      % mono
%      run everything below from "STAGE 1" on.
%  Without a file the script synthesizes a 60 s test signal
%  (8 power-law onsets + decay, DC, 50 Hz hum, 0.1 Hz drift, noise)
%  and scores itself against ground truth.
%
%  Verified (R2024a base MATLAB, this machine):
%    60 s / 8 onsets (A = 0.3..1.0, beta = 0.55..0.85, sigma = 0.02):
%      8/8 detected, 0 false alarms, |dT| <= 31 samples (0.7 ms),
%      1.7 s total (~35x realtime).
%    5-seed robustness: 39/40 detected, 0 false alarms in all 5,
%      max |dT| = 47 samples (1.1 ms).

clear; clc;

% Project root: hardcoded, as in the four replay drivers and demo.m
% (this environment's wrapper run() reports mfilename('fullFile') as a
% bare basename, so it cannot be derived portably from here).
root = 'D:/Source/hermes-dir/fdedge';
if ~isfile(fullfile(root, 'fdedge', 'config.m'))
    error('fdedge:sweep', 'project root not found at %s', root);
end
addpath(fullfile(root, 'fdedge'));
addpath(fullfile(root, 'fdedge', 'fd'));
addpath(fullfile(root, 'fdedge', 'est'));
addpath(fullfile(root, 'fdedge', 'calib'));
addpath(fullfile(root, 'fdedge', 'tracker'));
addpath(fullfile(root, 'fdedge', 'signal'));

%% ------------------------------------------------------------------ config
M      = 1024;      % window length in samples (23.2 ms @ 44.1 kHz).
                    %   Must match the calibration; rebuild calib per M (D7).
                    %   P4 guard: M <= 8192.
H      = M/4;       % hop (75% overlap). Edge must fall in the interior
                    %   [0.25, 0.75] zone; adjacent windows cover the rest.
Wbox   = 96;        % step-profile box half-width (samples)
dgap   = 32;        % gap between trial edge and boxes (samples): keeps the
                    %   boxes away from the response's near-edge region
thrH   = 0.2;       % stage-1 threshold on max step height
                    %   (measured noise: median 0.072, q99 0.126;
                    %    weakest verified onset ~0.37 -> 11 sigma)
dedupM = 0.3*M;     % final dedup: keep best candidate within this span
gateSNR = 3;        % stage-2 verification: stepSNR gate
                    %   (measured noise candidates: all <= 1.4)
gatePos = 0.0201;   % stage-2 position sanity: refine search width is 0.02
                    %   by construction, so this accepts every refine result
                    %   and is diagnostic; stepSNR is the real gate

%% ------------------------------------------------------------ test signal
% --- PLUG YOUR FILE HERE: [y, Fs] = audioread('clip.wav'); y = y(:,1); ---
if ~exist('y', 'var')
    rng(99);
    Fs  = 44100;
    N   = 60*Fs;
    n   = (0:N-1)';
    tau = 4410;   % 100 ms decay after the M-sample attack
    K   = 8;
    Ttrue = [3.2 7.8 12.4 18.9 25.1 33.6 44.2 55.7].*Fs;
    betat = [0.55 0.70 0.70 0.85 0.55 0.70 0.85 0.70];
    At    = [0.9 0.6 0.4 0.8 0.3 1.0 0.5 0.7];
    y = 0.3 + 0.01*sin(2*pi*50*n/Fs) + 0.05*sin(2*pi*0.1*n/Fs) ...
        + 0.02*randn(N,1);
    for k = 1:K
        d = n - Ttrue(k);
        % attack (power law over M samples) + exponential decay:
        % a REAL transient is bounded — the fdedge cusp is a LOCAL
        % attack model, so the synthetic must decay or sustain.
        y = y + At(k)*(max(d,0)/M).^betat(k) .* exp(-max(d - M, 0)/tau);
    end
    fprintf('synthetic file: %d samples (%.0f s), %d onsets\n', N, N/Fs, K);
end

%% ------------------------------------------------------------- stage 1
x     = (0:M)'/M;
nsamp = floor((N - M)/H);
hmax  = zeros(nsamp,1);
xloc  = zeros(nsamp,1);

t1 = tic;
for iw = 1:nsamp
    n0 = 1 + (iw-1)*H;
    w  = y(n0:n0+M);
    % required for real signals: remove the local DC level using only
    % pre-edge samples (left quarter)
    w  = w - median(w(1:round(0.25*(M+1))));
    fs = smoothGauss(w, 8);
    R  = caputoGL(x, fs, 0.70);          % one fixed detection order

    % step-height profile: H(p) = mean(R over [p+dgap, p+dgap+Wbox-1])
    %                      - mean(R over [p-dgap-Wbox+1, p-dgap])
    % via prefix sums (O(M) per window)
    cR = [0; cumsum(R)];
    ii = (dgap+Wbox) : (M+2-dgap-Wbox);
    Hv = (cR(ii+dgap+Wbox) - cR(ii+dgap) - cR(ii-dgap) ...
          + cR(ii-dgap-Wbox+1)) / Wbox;
    im = (ii >= 0.25*M) & (ii <= 0.75*M);   % interior zone only
    Hv2 = Hv; Hv2(~im) = -Inf;
    [hm, jk] = max(Hv2);
    hmax(iw) = hm;

    % locator: centroid of the positive profile around the peak.
    % For a true edge the profile is a trapezoid symmetric about the
    % edge -> centroid = edge (+/- a few samples). The slope argmax is
    % a BAD locator for beta > alpha edges (edge slope ~ 1 sigma).
    lo = max(1, jk - round(0.15*M));
    hi = min(numel(Hv), jk + round(0.15*M));
    Hc = max(Hv(lo:hi), 0);
    xc = ii(lo:hi); xc = xc(:);
    xloc(iw) = (sum(xc .* Hc) / max(sum(Hc), eps)) / M;
end
t1 = toc(t1);

% threshold + NMS in window units + dedup in sample units
iCand = find(hmax > thrH);
[sv, iord2] = sort(hmax(iCand), 'descend');
minGap = round(M/(2*H));
keptW = [];
for j = 1:numel(iCand)
    c = iCand(iord2(j));
    if isempty(keptW) || all(abs(c - keptW) >= minGap)
        keptW = [keptW; c];
    end
end
keptW = sort(keptW);
T0all = zeros(numel(keptW),1);
for k = 1:numel(keptW)
    iw = keptW(k);
    n0 = 1 + (iw-1)*H;
    T0all(k) = n0 - 1 + xloc(iw)*M;
end
[~, ord] = sort(T0all);
T0all = T0all(ord);
hord  = hmax(keptW(ord));
keep  = true(numel(T0all),1);
for a = 1:numel(T0all)
    for b = a+1:numel(T0all)
        if T0all(b) - T0all(a) < dedupM
            if hord(b) > hord(a), keep(a) = false;
            else                 , keep(b) = false; end
        end
    end
end
keptW = keptW(ord(keep));
T0all = T0all(keep);

fprintf('stage1 (%.2f s): hmax median %.3f, q99 %.3f, thr %.2f\n', ...
        t1, median(hmax), quantile(hmax,0.99), thrH);
fprintf('  -> %d raw windows, %d after NMS, %d after dedup\n', ...
        numel(iCand), numel(keptW), numel(T0all));

%% ------------------------------------------------------------- stage 2
cfg   = config('Profile','beta','M',M);
calib = gridSlopePeak(cfg);

t2 = tic;
npass = 0;
fprintf('\n%-3s %12s %9s %8s %8s %8s %8s  %s\n', ...
        'k','T0 (smp)','dTref','beta_h','A_h','stepSNR','hmax','verdict');
for k = 1:numel(keptW)
    iw = keptW(k);
    n0 = 1 + (iw-1)*H;
    T0 = T0all(k);
    w  = y(n0:n0+M);
    w  = w - median(w(1:round(0.25*(M+1))));
    xguess = (T0 - (n0-1))/M;

    trk = Tracker(cfg, calib);
    [st, ~] = trk.update(w, xguess);
    Th = (n0-1) + st.xraw*M;

    fs2 = smoothGauss(w, cfg.sGrid);
    R2  = caputoGL(x, fs2, st.braw);
    iE  = min(max(round(st.xraw*M), 1), M);
    stepSNR = mean(R2(iE+32:min(M+1, iE+128))) / std(R2(1:64));
    Ahat    = mean(R2(iE+52:min(M+1, iE+205)));

    ok = (stepSNR > gateSNR) && (abs(st.xraw - xguess) <= gatePos);
    if ok, npass = npass + 1; end

    if exist('Ttrue', 'var')
        im2 = find(abs(Ttrue - Th) < 0.6*M);
        if ~isempty(im2)
            dT = Th - Ttrue(im2(1));
            fprintf('%-3d %12.0f %+9.0f %8.4f %8.3f %8.2f %8.3f  %s\n', ...
                    k, T0, dT, st.braw, Ahat, stepSNR, hmax(iw), ...
                    tern(ok,'PASS','reject'));
        else
            fprintf('%-3d %12.0f %9s %8.4f %8.3f %8.2f %8.3f  %s nomatch\n', ...
                    k, T0, '-', st.braw, Ahat, stepSNR, hmax(iw), ...
                    tern(ok,'PASS','reject'));
        end
    else
        fprintf('%-3d %12.0f %9s %8.4f %8.3f %8.2f %8.3f  %s\n', ...
                k, T0, '-', st.braw, Ahat, stepSNR, hmax(iw), ...
                tern(ok,'PASS','reject'));
    end
end
t2 = toc(t2);
fprintf('\nstage2 (%.2f s): %d verified events\n', t2, npass);
fprintf('total %.2f s for %.0f s of signal (~%.0fx realtime)\n', ...
        t1+t2, N/Fs, (N/Fs)/(t1+t2));

function s = tern(c, a, b)
    if c, s = a; else, s = b; end
end
