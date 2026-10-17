function delta = scalarHalfCross(cfg)
% Scalar half-crossing bias: noise-free one-sided power-law edge at
% xref = 0.55, smoothed, L1 derivative at the true order beta; the bias
% is the half-crossing position error delta = xhat - 0.55. Shift
% invariance makes delta a constant (independent of the edge position),
% so one scalar calibrates the whole track.
% (Move of the Part B calibration block in FD_EdgeDetector_PoC/
% edge_fd_poc.m.)
% Required cfg fields: M, A, beta, sGrid.
M = cfg.M;
x = (0:M)'/M;
kw = (-3*cfg.sGrid:3*cfg.sGrid)';
kg = exp(-kw.^2/(2*cfg.sGrid^2)); kg = kg/sum(kg);
fe0 = cfg.A*max(x - 0.55, 0).^cfg.beta;
fs0 = conv(fe0, kg, 'same');
R0  = caputoL1(x, fs0, cfg.beta);
P0  = quantile(R0, 0.90);
delta = halfCross(R0, x, 0.5*P0) - 0.55;
end
