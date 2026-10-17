function e = slopePeak(R, xgrid, ell, xlo, xhi)
% Position of the maximum of the central-difference slope of R
% (width 2*ell), with parabolic sub-grid interpolation. For a power-law
% edge the D^alpha response has its steepest point exactly at the edge
% for ALL alpha in (0,1) (step for alpha=beta, power-law corner otherwise).
% Optional xlo/xhi restrict the search window (robustness vs noise spikes).
%
% Moved verbatim from FD_EdgeDetector_PoC/edge_fd_beta_track.m: slope_peak
% (including the max(abs(S)) fix from the PoC's history). Entry guards: R
% and xgrid coerced to columns (no-op for the PoC's column inputs).
% P3: the window restriction (absS = -Inf outside [xlo, xhi]) is what
% rejects the spurious zero-pad boundary edge at x = 1.
% NOTE: the PoC always feeds this a Gaussian-smoothed response; the
% unsmoothed GL response rings at a corner and the peak is biased ~ell
% grid (see test_slopePeak/unsmoothedRingingDocumented).
R = R(:);
xgrid = xgrid(:);
hh = xgrid(2) - xgrid(1);
S  = (R(1+2*ell:end) - R(1:end-2*ell)) / (2*ell*hh);
xs = xgrid(1+ell:end-ell);
absS = abs(S);
if nargin >= 5
    absS(xs < xlo | xs > xhi) = -Inf;
end
[sv, ik] = max(absS);
if ik <= 1 || ik >= numel(S)
    e = xs(ik);
    return;
end
d2 = S(ik-1) - 2*S(ik) + S(ik+1);
if abs(d2) < 1e-14
    e = xs(ik);
else
    e = xs(ik) + 0.5*(S(ik-1) - S(ik+1))/d2 * hh;
end
end
