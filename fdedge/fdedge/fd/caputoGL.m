function D = caputoGL(t, f, a)
% Caputo fractional derivative, Grunwald-Letnikov form via convolution,
% O(N log N).  D^a f(t_j) = h^(-a) * sum_k w(j-k+1) * (f(k) - f(1)),
% w(1) = 1, w(k+1) = w(k)*(k-1-a)/k.
% Moved verbatim from FD_EdgeDetector_PoC/edge_fd_beta_track.m: fdGL
% (rename + package prefix only; f = f(:) entry guard, D13/P4).
f = f(:);
N  = numel(t) - 1;
hh = t(2) - t(1);
w  = ones(N+1, 1);
for k = 1:N
    w(k+1) = w(k)*(k-1-a)/k;
end
Dfull = hh^(-a) * conv(f - f(1), w, 'full');
D = Dfull(1:numel(t));
end
