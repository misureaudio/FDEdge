%% ROBUST_STEP_SWEPT — co-gate novelty of a 0.5 step at 6 positions
%
% Robustness check for the spectral-novelty co-gate (est/spectralNovelty.m):
% the same 0.5 step injected at 6 different times on the PrimeTone
% background; the co-gate must pass every position (nov > 0.1).
% Stage 1 identical to examples/stress_multitone_gate.m.

root = 'D:/Source/hermes-dir/fdedge';
addpath(fullfile(root, 'fdedge'));
addpath(fullfile(root, 'fdedge', 'fd'));
addpath(fullfile(root, 'fdedge', 'est'));
addpath(fullfile(root, 'fdedge', 'calib'));
addpath(fullfile(root, 'fdedge', 'tracker'));
addpath(fullfile(root, 'fdedge', 'signal'));

flac = 'C:/Users/MATTIA/AppData/Local/hermes/attachments/PrimeToneNewPhi[Nch 1][Frac 03][LFL 043][00020-21600][Nfr 027]-[48000-24]-[1s]-[-0.25dB]-[049062]-[ 8.5905]-[RD].flac';
[y, Fs] = audioread(flac);
y = y(:,1);
N = numel(y);
t = ((0:N-1)')/Fs;

M = 1024; H = M/4; Wbox = 96; dgap = 32;
x = (0:M)'/M;
cfg = config('Profile','beta','M',M);
calib = gridSlopePeak(cfg);
nsamp = floor((N - M)/H);

tEdges = [0.15 0.30 0.45 0.60 0.75 0.90];
for k = 1:numel(tEdges)
    tEdge = tEdges(k);
    yA = y + 0.5*(t >= tEdge);
    sEdge = tEdge*Fs + 1;
    bestNov = 0; bestIw = 0; bestH = 0; bestIE = 0;
    for iw = 1:nsamp
        n0 = 1 + (iw-1)*H;
        if ~(sEdge >= n0 + 0.25*M && sEdge <= n0 + 0.75*M), continue; end
        w = yA(n0:n0+M);
        w = w - median(w(1:round(0.25*(M+1))));
        fs = smoothGauss(w, 8);
        R = caputoGL(x, fs, 0.70);
        cR = [0; cumsum(R)];
        ii = (dgap+Wbox):(M+2-dgap-Wbox);
        Hv = (cR(ii+dgap+Wbox) - cR(ii+dgap) - cR(ii-dgap) ...
              + cR(ii-dgap-Wbox+1)) / Wbox;
        im = (ii >= 0.25*M) & (ii <= 0.75*M);
        Hv2 = Hv; Hv2(~im) = -Inf;
        [hm, jk] = max(Hv2);
        lo = max(1, jk - round(0.15*M));
        hi = min(numel(Hv), jk + round(0.15*M));
        Hc = max(Hv(lo:hi), 0);
        xc = ii(lo:hi); xc = xc(:);
        xloc = (sum(xc .* Hc) / max(sum(Hc), eps)) / M;
        iE = min(max(round(xloc*M), 16), M-16);
        % co-gate with the same position scan as stress_multitone_gate.m
        nov = 0;
        for d = [-32 0 32 64]
            i = min(max(iE + d, 16), M+1-16);
            [nv, ~] = spectralNovelty(w, i, Fs);
            nov = max(nov, nv);
        end
        if nov > bestNov
            bestNov = nov; bestIw = iw; bestH = hm; bestIE = iE;
        end
    end
    tag = 'REJECT';
    if bestNov > 0.1, tag = 'PASS'; end
    fprintf('step t=%.2f s: best window %d (iE %d), hmax %.3f, nov %.4f (%s)\n', ...
            tEdge, bestIw, bestIE, bestH, bestNov, tag);
end
fprintf('ROBUST_DONE\n');
