function e = cuspMin(R, xgrid, xlo, xhi)
% Position of the minimum of R (parabolic sub-grid), search restricted to
% [xlo, xhi]. For a symmetric power-law edge the Riesz response
% (-Delta)^{a/2}|x-xe|^beta has a cusp-shaped local MINIMUM at the edge
% for all a in (0,1).
%
% Moved verbatim from FD_EdgeDetector_PoC/edge_fd_hybrid_track.m: argmin_cusp
% (identical copy in edge_fd_riesz_track.m). Entry guards: R and xgrid are
% coerced to columns (no-op for the PoC's column inputs).
% P3: the window restriction (Rw = Inf outside [xlo, xhi]) is what rejects
% the spurious zero-pad boundary edges; keep the window tight around the
% predicted position (the PoC uses +/-0.02 around xguess).
R = R(:);
xgrid = xgrid(:);
hh = xgrid(2) - xgrid(1);
Rw = R;
Rw(xgrid < xlo | xgrid > xhi) = Inf;
[~, ik] = min(Rw);
if ik <= 1 || ik >= numel(R)
    e = xgrid(ik);
    return;
end
y0 = Rw(ik-1); y1 = Rw(ik); y2 = Rw(ik+1);
dd = y0 - 2*y1 + y2;
if ~isfinite(dd) || abs(dd) < 1e-14
    e = xgrid(ik);
else
    e = xgrid(ik) + (y0 - y2)/(2*dd)*hh;
end
end
