function [bhat, iaStar] = exponentBank(c1, alphas)
% Zero of the c1-vs-alpha line, clamped to the bank; iaStar = the bank
% member closest to bhat. Theory: c1(a) = g*(beta - a) with g > 0
% depending only on the window geometry, so the zero of the line is beta
% (invariant to the window placement; the raw intercept alone carries a
% systematic bias proportional to the xguess/xe offset).
% Moved verbatim from FD_EdgeDetector_PoC/edge_fd_beta_track.m /
% edge_fd_hybrid_track.m (est_exponent tail):
%   p = polyfit(alphas, c1, 1);
%   if abs(p(1)) > 1e-12
%       bhat = min(max(-p(2)/p(1), alphas(1)), alphas(end));
%   else
%       bhat = alphas(round(numel(alphas)/2));
%   end
%   [~, ia_star] = min(abs(alphas - bhat));
% Zero slope -> the bank's center member alphas(round(numel(alphas)/2))
% (1-based: for the 7-member PoC bank that is 0.70, not 0.75).
c1 = c1(:);
alphas = alphas(:);
p = polyfit(alphas, c1, 1);
if abs(p(1)) > 1e-12
    bhat = min(max(-p(2)/p(1), alphas(1)), alphas(end));
else
    bhat = alphas(round(numel(alphas)/2));
end
[~, iaStar] = min(abs(alphas - bhat));
end
