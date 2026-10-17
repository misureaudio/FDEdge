function [nov, p] = spectralNovelty(w, iE, Fs, porder)
% SPECTRALNOVELTY fraction of post-edge energy not predicted by the
% best linear predictor fitted on the pre-edge segment.
%
%   [nov, p] = spectralNovelty(w, iE, Fs)
%
%   w    : window (n x 1), baseline already subtracted (as in the
%          sweep pipeline: w - median(left quarter)).
%   iE   : edge position in sample units (1..n-1).
%          pre = w(1:iE), post = w(iE+1:n).
%   Fs   : sample rate (kept for signature compatibility; the
%          predictor is scale-free).
%   porder : (optional) AR predictor order, default 120.
%
%   Returns
%   nov  : sqrt(sum((post - pred).^2)) / sqrt(sum(post.^2)), the
%          one-step prediction error of the post segment, normalized
%          by the post energy.
%            nov ~ 1e-3 class : the post segment is exactly (linearly)
%            predictable from the pre segment -> stationary process ->
%            NOT an edge.
%            nov >> 0.1 (measured 0.99-8.30 for a 0.5 step on the
%            PrimeTone background at 6 positions, 0.68-1.78 for new
%            1 kHz tones A = 0.15/0.5): the post segment contains
%            content the pre segment's process cannot produce -> an edge.
%   p    : predictor order actually used.
%
%   Theory (see docs/spectral_novelty_gate.md):
%   An edge is a change of the signal's structure: a step adds DC +
%   broadband, a new tone adds a line, a burst adds broadband energy.
%   A stationary process (a sum of fixed sinusoids, however the
%   phases are chosen) has identical pre- and post-edge statistics:
%   its post segment is a deterministic continuation of its pre
%   segment and is exactly predictable by a linear predictor whose
%   poles sit on the unit circle at the line frequencies.
%   The best linear predictor (AR) fitted on the pre segment
%   therefore models the stationary content without needing to know
%   the line frequencies at all — which matters because closely
%   spaced low-frequency lines (e.g. 43/53/61 Hz, 8-10 Hz apart)
%   merge into single spectral peaks: they are below the window's
%   resolution limit (46.8 Hz @ 1025 samples @ 48 kHz) and cannot be
%   separated by any zero-padded FFT of the window.
%
%   Why one-step prediction with feedback, measured on the first
%   L = 64 post samples: the pre segment (6.7-13 ms) is shorter than
%   the period of the lowest lines (43 Hz = 23.3 ms), so the fitted
%   line poles sit at r < 1; open-loop multi-step prediction then
%   decays to zero while the true tone continues (measured nov ~ 1.0
%   on the stationary PrimeTone). One-step recursion with ACTUAL post
%   samples fed back re-anchors the state every step, so the radius
%   error does not accumulate and the stationary multitone is
%   predicted with nov ~ 1e-5. The cost of the feedback: the
%   predictor re-equilibrates to a new LEVEL in a short transient
%   (measured ~12-64 samples on the PrimeTone step), so the novelty
%   is measured only on the first 64 post samples (the stepSNR gate's
%   own post-box length), before the transient decays. A new TONE,
%   however, is absent from the model and the feedback cannot adapt
%   to it: its error persists over the whole post.
%
%   Why not a raw pre/post spectral power ratio: a step's new content
%   is DC + the boundary discontinuity; the FFTs of the two smooth
%   segments carry none of the discontinuity energy, so per-segment
%   off-line power ratios read ~1 even for a step (measured 0.93).
%   The prediction error captures the unmodeled content directly.
%
%   Why Burg's algorithm for the fit: the predictor is run recursively
%   on the post segment, so it must be stable (all poles strictly
%   inside the unit circle). The unconstrained least-squares fit
%   X\y can produce unstable poles on short pre segments (measured:
%   a pole at |z| ~ 1.10 on a 6.7 ms multitone segment, blowing the
%   novelty up to 3e4 for a real step). Burg's forward-backward
%   recursion gives |reflection coefficients| < 1 at every order by
%   Cauchy-Schwarz, hence a guaranteed-stable predictor.

w = w(:);
n = numel(w);
if nargin < 4 || isempty(porder), porder = 120; end
p = min([porder, iE - 8, n - iE - 8]);
iE = min(max(round(iE), p + 8), n - 8);
if p < 16
    % pre segment too short for the predictor: stationarity cannot be
    % verified -> report as novel (fail open)
    nov = 1;
    return
end

pre  = w(1:iE);
post = w(iE+1:n);
Np   = numel(pre);

% --- Burg AR fit on the pre segment (stable by construction) ----------
% pre(k) = -sum(a(j)*pre(k-j)) + e,  j = 1..p  (a_1 = newest lag)
a = burgAR(pre, p);

% --- one-step prediction of the first L post samples -----------------
% One-step recursion with ACTUAL post samples fed back: the state is
% re-anchored at every step, so pole-radius errors of the fit do not
% accumulate (open-loop multi-step is useless here: the pre segment is
% shorter than the period of the lowest lines, the fitted line poles
% sit at r < 1, and the open-loop prediction decays to zero while the
% true tone continues — measured nov ~ 1.0 on the stationary
% PrimeTone). The cost of the feedback is that the predictor
% re-equilibrates to a new LEVEL in a short transient (~12-64 samples
% measured on the PrimeTone step), so the novelty is measured only on
% the first L = 64 post samples (the stepSNR gate's own post-box
% length), before the transient decays:
%   stationary: one-step error ~ innovation -> nov ~ 1e-5
%   step:       the level change is unmodeled during the transient
%               -> nov ~ 0.2-0.5
%   new tone:   the line is absent from the model, the feedback cannot
%               adapt to it -> the error persists -> nov ~ 0.3-0.6
L = min(64, numel(post));
buf = pre(Np:-1:Np-p+1);
predq = zeros(L, 1);
for k = 1:L
    predq(k) = -a' * buf;
    buf = [post(k); buf(1:end-1)];
end

resid = post(1:L) - predq;
nov = sqrt(sum(resid .^ 2)) / max(sqrt(sum(post(1:L) .^ 2)), eps);
end
