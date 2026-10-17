function c1 = loglogSlope(y, dvec)
% Log-log slope of y over the distance vector dvec, with the PoC's 1e-6
% floor on y (the polyfit pattern used by all four PoCs):
%   p = polyfit(log(dvec), log(max(y, 1e-6)), 1); c1 = p(1);
% y: response samples (vector); dvec: corresponding distances (vector, > 0).
% The 1e-6 floor keeps all-zero / near-zero responses finite (slope ~ 0).
% Moved verbatim from FD_EdgeDetector_PoC (the polyfit pattern).
% Entry guards: y and dvec coerced to columns.
y = y(:);
dvec = dvec(:);
pf = polyfit(log(dvec), log(max(y, 1e-6)), 1);
c1 = pf(1);
end
