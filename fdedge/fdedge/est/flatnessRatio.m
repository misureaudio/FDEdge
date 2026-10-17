function [ratio, c1, platMean, iBest] = flatnessRatio(R, alphas, N, w, iBeta)
% Part A flatness metrics of FD_EdgeDetector_PoC/edge_fd_poc.m, moved
% verbatim (the ratio/plateau block, lines ~45-58):
%   ratio(j)   = mean(R(k_hi-w+1:k_hi, j)) / mean(R(k_lo:k_lo+w-1, j))
%                (endpoint ratio R(0.95)/R(0.05); = 1 only at the matched
%                order alpha = beta)
%   c1(j)      = log-log slope of R(k_lo:k_hi, j) vs the index
%                (estimates beta - alpha)
%   platMean   = mean(R(k_lo:k_hi, iBeta))  (plateau at the matched order;
%                = A*gamma(beta+1) for f = A*t^beta)
%   iBest      = argmin_j |ratio(j) - 1|  (best-matching order)
% R: (N+1) x numel(alphas) L1 response bank on t = (0:N)'; alphas: column.
% NOTE: the plan's API listed [ratio, platMean, iBest]; c1 is returned as
% a 4th output (computed in the same verbatim loop; needed by the Part A
% replay driver).
k_lo = round(0.05*N)+1;
k_hi = round(0.95*N)+1;
ratio = zeros(size(alphas));
c1    = zeros(size(alphas));
for j = 1:numel(alphas)
    ratio(j) = mean(R(k_hi-w+1:k_hi, j)) / mean(R(k_lo:k_lo+w-1, j));
    yfull = R(k_lo:k_hi, j);
    p = polyfit(log(k_lo:k_hi), log(max(yfull, 1e-6)), 1);
    c1(j) = p(1);
end
platMean = mean(R(k_lo:k_hi, iBeta));
[~, iBest] = min(abs(ratio-1));
end
