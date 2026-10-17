function R = rieszLaplacian(f, a, dx)
% Windowed-FFT fractional Laplacian (-Delta)^{a/2}; scalar or bank orders.
% k MUST be a column (implicit-expansion trap, P4). Moved verbatim from
% FD_EdgeDetector_PoC/edge_fd_riesz_track.m: riesz_bank (loop body
% unchanged; scalar/bank dispatch + f = f(:) + M guard added, D13/P4).
% dx: grid spacing, optional. Default 1/M = the PoC period-1 convention
% (k = 2*pi*k_idx, byte-identical to riesz_bank). Pass the physical
% spacing to use the physical wavenumbers 2*pi*k_idx/(M*dx).
% ACCURACY NOTE: this computes the PERIOD-L operator (L = M*dx). For
% 0 < a < 1 the |k|^a multiplier has a cusp at k=0, so the result
% approaches the continuum fractional Laplacian only algebraically,
% ~L^{-1-a}, as L grows (measured: error ratio 2^1.7 per doubling of L
% at a=0.7). Refining dx at fixed L does NOT reduce this error — grow
% the period (and M with it).
% k MUST be a column (implicit-expansion trap, P4). kidx is the signed FFT
% frequency order [0:floor(M/2), -floor((M-1)/2):-1], built for any M
% (even or odd); for even M it is identical to the PoC riesz_bank vector
% [0:M/2, -(M/2-1):-1], so PoC byte-reproducibility is preserved.
f = f(:);
M = numel(f);
if M > 8192
    error('fdedge:bigM', 'rieszLaplacian: M = %d exceeds the 8192-point sandbox-safe limit', M);
end
if nargin < 3 || isempty(dx)
    dx = 1/M;
end
kidx = (0:M-1) - M*((0:M-1) > M/2);
k = (kidx * (2*pi/(M*dx)))';
F = fft(f);
if isscalar(a)
    R = real(ifft(F .* abs(k).^a));
else
    R = zeros(M, numel(a));
    for ia = 1:numel(a)
        R(:,ia) = real(ifft(F .* abs(k).^a(ia)));
    end
end
end
