function fs = smoothGauss(f, sGrid)
% Gaussian pre-smoothing: normalized 6-sigma kernel, conv(..., 'same').
% Moved verbatim from the kw/kg + conv pattern in all four PoC scripts
% (edge_fd_poc.m:25-26 et al.; rename + package prefix only; f = f(:)
% entry guard, D13/P4).
% P3: zero-padded 'same' smoothing creates a spurious boundary edge (the
% smoothed signal is forced toward zero at the domain ends) - all
% extremum searches must be window-restricted (see est.* docstrings).
f = f(:);
kw = (-3*sGrid:3*sGrid)';
kg = exp(-kw.^2/(2*sGrid^2)); kg = kg/sum(kg);
fs = conv(f, kg, 'same');
end
