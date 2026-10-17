function a = burgAR(x, p)
% BURGAR  Burg's AR(p) estimate (forward-backward Levinson recursion).
%
%   a = burgAR(x, p) returns the AR coefficients a (p x 1, newest lag
%   first) of the model
%
%       x(k) = -sum(a(j)*x(k-j)) + e(k),  j = 1..p,
%
%   fitted on the raw data (no mean removal: a DC level is a
%   legitimate pole of the process and must be predicted, not
%   subtracted away).
%
%   Recursion (exact identities, with f_m the forward and b_m the
%   backward prediction error of order m, valid for n = m+1..N):
%
%       f_m(n) = f_{m-1}(n) + k_m b_{m-1}(n-1)
%       b_m(n) = b_{m-1}(n-1) + k_m f_{m-1}(n)
%
%   together with the Levinson coefficient update
%   a_m(i) = a_{m-1}(i) + k_m a_{m-1}(m-i), a_m(m) = k_m.
%
%   k_m is chosen to minimize the combined forward + backward error
%   energy of order m (each cross term appears in both brackets,
%   hence the factor 2):
%
%       k_m = -2 sum f_{m-1}(n) b_{m-1}(n-1)
%             / (sum f_{m-1}(n)^2 + sum b_{m-1}(n-1)^2),  n = m+1..N
%
%   so |k_m| < 1 at every step by Cauchy-Schwarz + AM-GM (strict
%   unless the two error sequences are proportional) -> all poles of
%   the resulting predictor are strictly inside the unit circle ->
%   the one-step recursion is stable for ANY input. That stability
%   guarantee is what the spectral-novelty co-gate needs on short
%   pre segments, where the unconstrained least-squares fit X\y can
%   return unstable poles (measured |z| ~ 1.10 on a 6.7 ms multitone
%   segment; the recursive novelty then blew up to 3e4).
%
%   Requires numel(x) >= p + 2.

    x = x(:);
    N = numel(x);
    assert(N >= p + 2, 'burgAR: need numel(x) >= p + 2, got %d and p = %d', N, p);
    a = zeros(p, 1);
    den = sum(x(1:N-1).^2) + sum(x(2:N).^2);
    k = -2*sum(x(1:N-1).*x(2:N)) / max(den, eps);
    a(1) = k;
    f = x(2:N) + k*x(1:N-1);      % f_1(n), n = 2..N
    b = x(1:N-1) + k*x(2:N);      % b_1(n), n = 2..N
    for m = 2:p
        num  = sum(f(2:end).*b(1:end-1));
        den  = sum(f(2:end).^2) + sum(b(1:end-1).^2);
        k = -2*num / max(den, eps);
        aold = a(1:m-1);
        a(1:m-1) = aold + k*aold(m-1:-1:1);
        a(m) = k;
        fn = f(2:end) + k*b(1:end-1);   % f_m(n), n = m+1..N
        bn = b(1:end-1) + k*f(2:end);   % b_m(n), n = m+1..N
        f = fn; b = bn;
    end
end
