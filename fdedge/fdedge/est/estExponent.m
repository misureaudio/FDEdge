function [bhat, iaStar] = estExponent(fs, alphas, xref, x, w0, w1)
% Memory-subtracted one-sided Caputo exponent at position xref.
% g = fs with the left side (x <= xref) zeroed -> D^a g is the right-side
% cusp only (the left-side memory term is removed exactly, with no extra
% noise: g is exactly zero on the left). c1(a) = log-log slope of D^a g
% over [xref+w0, xref+w1]; bhat = zero of the c1-vs-a line, clamped to
% the bank; iaStar = the bank member closest to bhat.
% (Move of the local function in FD_EdgeDetector_PoC/
% edge_fd_hybrid_track.m, composed from the M4 primitives
% memorySubtract / caputoGL / loglogSlope / exponentBank - bit-identical
% to the PoC's est_exponent.)
M  = numel(x) - 1;
g  = memorySubtract(fs, xref, M);
ka = min(M+1, max(2, round((xref + w0)*M) + 1));
kb = min(M+1, max(ka+2, round((xref + w1)*M) + 1));
dvec = x(ka:kb) - xref;
c1 = zeros(size(alphas));
for ia = 1:numel(alphas)
    D = caputoGL(x, g, alphas(ia));
    c1(ia) = loglogSlope(D(ka:kb), dvec);
end
[bhat, iaStar] = exponentBank(c1, alphas);
end
