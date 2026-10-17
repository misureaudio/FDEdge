function D = caputoL1(t, f, a)
% L1 (product-integration) Caputo fractional derivative, uniform grid.
%  D^a f(t_j) = sum_k b_{j,k} * (f(k+1)-f(k))/h
%  b_{j,k} = ((t_j - t_k)^(1-a) - (t_j - t_{k+1})^(1-a))/gamma(2-a)
% Moved verbatim from FD_EdgeDetector_PoC/edge_fd_poc.m: fdL1
% (rename + package prefix only; f = f(:) entry guard, D13/P4).
% P1 guard (the one permitted deviation): the L1 scheme degenerates at
% alpha = 1 (silent zeros) - reject alpha >= 1 explicitly.
if a >= 1 - eps
    error('fdedge:badAlpha', 'caputoL1 is defined for alpha in (0,1); got %.15g', a);
end
f = f(:);
N  = numel(t) - 1;
hh = t(2) - t(1);
d  = diff(f);
D  = zeros(size(t));
gam = gamma(2 - a);
for jj = 1:N
    b = ((t(jj+1) - t(1:jj)).^(1-a) - (t(jj+1) - t(2:jj+1)).^(1-a)) / gam;
    D(jj+1) = sum(b .* d(1:jj)) / hh;
end
end
