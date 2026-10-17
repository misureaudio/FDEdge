function e = halfCross(R, xgrid, level)
% First half-threshold crossing of the step-like response R, subgrid by
% linear interpolation; requires a run of 20 consecutive above-thr samples
% (rejects isolated noise spikes).
%
% Moved verbatim from FD_EdgeDetector_PoC/edge_fd_poc.m: cross_est.
% Entry guards: R and xgrid coerced to columns (no-op for the PoC's
% column inputs).
% Shift invariance: on a uniform grid the bias (halfCross - xe) is EXACTLY
% constant for edges sitting on grid points (exact integer translation,
% verified to 1e-16 in the test) and varies sub-grid for off-grid edges
% (their fractional cell position differs); the pocB single-reference
% calibration relies on this property.
R = R(:);
xgrid = xgrid(:);
above = R > level;
run = false(size(above));
cnt = 0;
for k = numel(above):-1:1
    if above(k)
        cnt = cnt + 1;
    else
        cnt = 0;
    end
    if cnt >= 20
        run(k) = true;
    end
end
istar = find(run, 1, 'first');
if isempty(istar)
    istar = find(above, 1, 'first');
end
if istar <= 1
    e = xgrid(1);
else
    e = xgrid(istar-1) + (xgrid(istar) - xgrid(istar-1)) * ...
        (level - R(istar-1)) / (R(istar) - R(istar-1));
end
end
