%% STRESS_MULTITONE_GATE — spectral-novelty co-gate on the PrimeTone file
%
% Stress test for the spectral-novelty co-gate (est/spectralNovelty.m)
% against the verified false-positive mode of the v1 gate: the stationary
% 27-tone PrimeTone multitone (docs/stationary_multitone_gate.md).
%
% Design (stage 1 unchanged from examples/sweep_demo.m):
%   stage 1   per-window (M=1024, H=M/4) Caputo-GL at alpha=0.70 +
%             step-height profile -> hmax, centroid edge position
%   co-gate   nov = coGate(w, iE, Fs): max spectralNovelty over the
%             position scan {iE-32, iE, iE+32, iE+64} (stage-1 centroid
%             jitter is measured up to ~90 samples early on this
%             background); a candidate is REJECTED if nov <= 0.1
%   old gate  (comparison only) stage-2 Tracker + stepSNR > 3 + position
%             sanity <= 0.0201
%
% Self-scoring (expect):
%   [x] co-gate false alarms on the stationary file: 0
%       (old gate measured 2026-10-17: 24 stepSNR passes, 20 with the
%        position co-gate)
%   [x] injected 0.5 step on the same background: co-gate PASS (nov > 0.1)
%   [x] injected new 1000 Hz tone (A=0.5) at 0.6 s: co-gate PASS
%       (the old gate is BLIND to it: stepSNR 0.48 — the derivative of
%        a 1 kHz tone oscillates ~20x inside the 96-sample box, box
%        mean ~ 0; this is the gap the co-gate fills)
%   [i] injected weak tone (A=0.15, -15 dB): old gate rejects
%       (stepSNR 1.08), co-gate PASSes (nov 0.68) — the co-gate is more
%       sensitive to weak tonal edges; in an AND combination with the
%       old gate the candidate is still rejected (sensitivity floor)

root = 'D:/Source/hermes-dir/fdedge';
addpath(fullfile(root, 'fdedge'));
addpath(fullfile(root, 'fdedge', 'fd'));
addpath(fullfile(root, 'fdedge', 'est'));
addpath(fullfile(root, 'fdedge', 'calib'));
addpath(fullfile(root, 'fdedge', 'tracker'));
addpath(fullfile(root, 'fdedge', 'signal'));

flac = 'C:/Users/MATTIA/AppData/Local/hermes/attachments/PrimeToneNewPhi[Nch 1][Frac 03][LFL 043][00020-21600][Nfr 027]-[48000-24]-[1s]-[-0.25dB]-[049062]-[ 8.5905]-[RD].flac';
assert(isfile(flac), 'PrimeTone FLAC not found at %s', flac);
[y, Fs] = audioread(flac);
y = y(:,1);
N = numel(y);
t = ((0:N-1)')/Fs;
fprintf('loaded %d samples @ %d Hz (%.3f s), peak %.4f, crest %.2f dB\n', ...
        N, Fs, N/Fs, max(abs(y)), 20*log10(max(abs(y))/rms(y)));

%% ------------------------------------------------------------ config
M      = 1024;      % window (matches the calibration)
H      = M/4;       % hop (75% overlap)
Wbox   = 96;        % step-profile box half-width
dgap   = 32;        % box-to-edge gap
thrH   = 0.2;     % stage-1 threshold
gateSNR = 3;      % old gate: stepSNR
gatePos = 0.0201;   % old gate: position sanity
gateNov = 0.1;      % co-gate: novelty floor

x     = (0:M)'/M;
cfg   = config('Profile','beta','M',M);
calib = gridSlopePeak(cfg);

%% ------------------------------------------- run 1: pure stationary file
t0 = tic;
[hmax, xloc, novs, snrs, posok] = runSweep(y, Fs, M, H, Wbox, dgap, ...
                                           x, cfg, calib, gatePos);
t1 = toc(t0);

fprintf('\nRUN 1 - pure stationary PrimeTone (no edge anywhere)\n');
fprintf('  stage1: %d windows, hmax median %.3f, q99 %.3f, > %.2f: %d\n', ...
        numel(hmax), median(hmax), quantile(hmax,0.99), thrH, sum(hmax > thrH));
fprintf('  co-gate nov: min %.3e, median %.3e, q99 %.3e, max %.3e\n', ...
        min(novs), median(novs), quantile(novs,0.99), max(novs));
nFA = sum(novs > gateNov);
fprintf('  co-gate (nov > %.2f): %d / %d windows pass -> %d FALSE ALARMS\n', ...
        gateNov, nFA, numel(novs), nFA);
nOld = sum(snrs > gateSNR);
nOldP = sum((snrs > gateSNR) & posok);
fprintf('  old gate (stepSNR > %d): %d / %d, with position co-gate: %d\n', ...
        gateSNR, nOld, numel(snrs), nOldP);
[sn, so] = sort(novs, 'descend');
fprintf('  top-10 nov windows: ');
for k = 1:10, fprintf('%.2e ', sn(k)); end
fprintf('\n');
ok1 = (nFA == 0);

%% ----------------------------- run 2: injected 0.5 step at t = 0.3 s
tEdge = 0.3;
yA = y + 0.5*(t >= tEdge);
[hmaxA, xlocA, novsA, snrsA, posokA] = runSweep(yA, Fs, M, H, Wbox, dgap, ...
                                                x, cfg, calib, gatePos);
sEdge = tEdge*Fs + 1;
% windows whose interior zone contains the edge
imW = false(numel(hmaxA),1);
for iw = 1:numel(hmaxA)
    n0 = 1 + (iw-1)*H;
    imW(iw) = (sEdge >= n0 + 0.25*M) && (sEdge <= n0 + 0.75*M);
end
idxA = find(imW);
[kA, jA] = max(novsA(idxA));
iwA = idxA(jA);
fprintf('\nRUN 2 - injected 0.5 step at t = %.2f s (same background)\n', tEdge);
fprintf('  best window %d: hmax %.3f, nov %.4f, stepSNR %.2f\n', ...
        iwA, hmaxA(iwA), kA, snrsA(iwA));
ok2 = (kA > gateNov);

%% ------------------------- run 3: injected new 1000 Hz tone at t = 0.6 s
% strong tone (A = 0.5). The OLD gate is blind to this edge: its
% stepSNR boxes the derivative profile, and the derivative of a 1 kHz
% tone oscillates ~20 times inside the 96-sample box, so the box MEAN
% ~ 0 (measured stepSNR = 0.48 << 3). The co-gate must pass it — this
% is the gap the co-gate exists to fill.
tEdge3 = 0.6;
yB = y + 0.5*cos(2*pi*1000*t).*(t >= tEdge3);
[hmaxB, xlocB, novsB, snrsB, posokB] = runSweep(yB, Fs, M, H, Wbox, dgap, ...
                                                 x, cfg, calib, gatePos);
sEdge3 = tEdge3*Fs + 1;
imW3 = false(numel(hmaxB),1);
for iw = 1:numel(hmaxB)
    n0 = 1 + (iw-1)*H;
    imW3(iw) = (sEdge3 >= n0 + 0.25*M) && (sEdge3 <= n0 + 0.75*M);
end
idxB = find(imW3);
[kB, jB] = max(novsB(idxB));
iwB = idxB(jB);
fprintf('\nRUN 3 - injected new 1000 Hz tone (A=0.5) at t = %.2f s\n', tEdge3);
fprintf('  best window %d: hmax %.3f, nov %.4f, stepSNR %.2f (old gate blind to new HF tones)\n', ...
        iwB, hmaxB(iwB), kB, snrsB(iwB));
ok3 = (kB > gateNov);

%% ------------------- run 4 (informational): weak tone A = 0.15
% a tone 15 dB below the background: the old gate rejects it (stepSNR
% 1.08 < 3); the co-gate PASSes it (nov 0.68 — the tone is a new line
% absent from the pre model, and its energy fraction over the 64-sample
% box exceeds the floor). The co-gate is thus MORE sensitive to weak
% tonal edges than the old gate; in an AND combination with the old
% gate the candidate is still rejected (sensitivity floor).
yC = y + 0.15*cos(2*pi*1000*t).*(t >= tEdge3);
[~, ~, novsC, snrsC, ~] = runSweep(yC, Fs, M, H, Wbox, dgap, ...
                                   x, cfg, calib, gatePos);
idxC = find(imW3);
[kC, jC] = max(novsC(idxC));
iwC = idxC(jC);
fprintf('\nRUN 4 - weak tone (A=0.15, -15 dB below background) at t = %.2f s\n', tEdge3);
fprintf('  best window %d: nov %.4f, stepSNR %.2f -> old gate %s, co-gate %s\n', ...
        iwC, kC, snrsC(iwC), ...
        tern(snrsC(iwC) > gateSNR, 'PASS', 'reject'), ...
        tern(kC > gateNov, 'PASS', 'reject'));

%% ---------------------------------------------------------------- score
fprintf('\nSCORE: run1 %s, run2 %s, run3 %s\n', ...
        tern(ok1, 'PASS', 'FAIL'), tern(ok2, 'PASS', 'FAIL'), ...
        tern(ok3, 'PASS', 'FAIL'));
fprintf('STRESS_MULTITONE_DONE (%.1f s for 4 sweeps of 1 s)\n', t1*4);

function s = tern(c, a, b)
    if c, s = a; else, s = b; end
end

function [hmax, xloc, novs, snrs, posok] = runSweep(y, Fs, M, H, Wbox, dgap, ...
                                                    x, cfg, calib, gatePos)
    N = numel(y);
    nsamp = floor((N - M)/H);
    hmax  = zeros(nsamp,1);
    xloc  = zeros(nsamp,1);
    novs  = zeros(nsamp,1);
    snrs  = zeros(nsamp,1);
    posok = false(nsamp,1);
    for iw = 1:nsamp
        n0 = 1 + (iw-1)*H;
        w  = y(n0:n0+M);
        w  = w - median(w(1:round(0.25*(M+1))));
        fs = smoothGauss(w, 8);
        R  = caputoGL(x, fs, 0.70);
        cR = [0; cumsum(R)];
        ii = (dgap+Wbox) : (M+2-dgap-Wbox);
        Hv = (cR(ii+dgap+Wbox) - cR(ii+dgap) - cR(ii-dgap) ...
              + cR(ii-dgap-Wbox+1)) / Wbox;
        im = (ii >= 0.25*M) & (ii <= 0.75*M);
        Hv2 = Hv; Hv2(~im) = -Inf;
        [hm, jk] = max(Hv2);
        hmax(iw) = hm;
        lo = max(1, jk - round(0.15*M));
        hi = min(numel(Hv), jk + round(0.15*M));
        Hc = max(Hv(lo:hi), 0);
        xc = ii(lo:hi); xc = xc(:);
        xloc(iw) = (sum(xc .* Hc) / max(sum(Hc), eps)) / M;
        % --- co-gate at the stage-1 edge position ---
        iE = min(max(round(xloc(iw)*M), 16), M-16);
        novs(iw) = coGate(w, iE, Fs);
        % --- old gate (stage 2, comparison only) ---
        trk = Tracker(cfg, calib);
        [st, ~] = trk.update(w, xloc(iw));
        fs2 = smoothGauss(w, cfg.sGrid);
        R2  = caputoGL(x, fs2, st.braw);
        iE2 = min(max(round(st.xraw*M), 1), M);
        snrs(iw)  = mean(R2(iE2+32:min(M+1, iE2+128))) / max(std(R2(1:64)), eps);
        posok(iw) = abs(st.xraw - xloc(iw)) <= gatePos;
    end
end

function nov = coGate(w, iE, Fs)
%COGATE spectral-novelty co-gate: max novelty over a small position
%scan around the stage-1 edge position.
%
%Why the scan: the stage-1 centroid is not exact on the multitone
%background (measured jitter: up to ~90 samples EARLY, ~30 late,
%depending on the local spectral mix). The novelty box is the first 64
%post samples; if the centroid is early by more than 64, the box falls
%entirely before the edge and a real edge reads nov ~ 1e-5 (measured:
%0.5 step at t = 0.90 s, centroid 88 early -> nov 0.0000). Scanning
%{iE-32, iE, iE+32, iE+64} covers centroid errors up to ~128 samples
%early and ~32 late, at the cost of 4 novelty evaluations per window.
best = 0;
for d = [-32 0 32 64]
    i  = min(max(iE + d, 16), numel(w) - 16);
    [nv, ~] = spectralNovelty(w, i, Fs);
    best = max(best, nv);
end
nov = best;
end
