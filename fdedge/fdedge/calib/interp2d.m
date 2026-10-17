function d = interp2d(G, axisRows, axisCols, qRow, qCol)
% Bilinear interpolation of a 2-D grid G at (qRow, qCol).
% axisRows = G's row axis (e.g. xes), axisCols = G's column axis (e.g.
% betas); both sorted ascending. Out-of-range queries are clipped to the
% range edges (no extrapolation).
% Moved verbatim from FD_EdgeDetector_PoC/edge_fd_riesz_track.m /
% edge_fd_hybrid_track.m: idx_w + the bilinear combination
%   dbr = w2*(w1*dbeta(i1,j1) + (1-w1)*dbeta(i1,j2)) + ...
%       (1-w2)*(w1*dbeta(i2,j1) + (1-w1)*dbeta(i2,j2));
%
% WEIGHT CONVENTION (verbatim PoC behavior): the combination
% w*G(i1) + (1-w)*G(i2) with w = (val - v(i1))/(v(i2) - v(i1)) evaluates
% the bilinear surface at the query point REFLECTED through the cell
% center. Exact at cell centers and at the range edges (i1 == i2); at an
% exact interior node (i,j) it returns the upper-right neighbor's value
% G(i+1,j+1). The calibration grids are smooth, so the deviation from
% standard bilinear is O(cell) and negligible; the M8 golden replay
% requires this verbatim behavior - do NOT "fix" the weights.
[i1, i2, wr] = idxw(axisRows, qRow);
[j1, j2, wc] = idxw(axisCols, qCol);
d = wc*(wr*G(i1,j1) + (1-wr)*G(i1,j2)) + (1-wc)*(wr*G(i2,j1) + (1-wr)*G(i2,j2));
end

function [i1, i2, w] = idxw(v, val)
% Bilinear index/weight of val in sorted vector v (clipped to range).
val = min(max(val, v(1)), v(end));
i = find(v <= val, 1, 'last');
if i >= numel(v)
    i1 = numel(v); i2 = numel(v); w = 0;
else
    i1 = i; i2 = i+1; w = (val - v(i1))/(v(i2) - v(i1));
end
end
