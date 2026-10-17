function g = memorySubtract(fs, xref, M)
% Memory subtraction (P6): g = fs with the left side (x <= xref) zeroed,
% so D^a g is the right-side cusp only - the left-side memory term of the
% Caputo derivative is removed EXACTLY, with no added noise (g is exactly
% zero on the left).
% Moved verbatim from FD_EdgeDetector_PoC/edge_fd_hybrid_track.m:
% est_exponent:  it = min(max(round(xref*M) + 1, 2), M); g = fs; g(1:it) = 0;
% M = numel(x) - 1 on x = (0:M)'/M (grid-interval count).
% Identity: D^a f - D^a f_left == D^a g (GL/L1 linearity; measured
% 2.8e-14, see test_memorySubtract).
g = fs(:);
it = min(max(round(xref*M) + 1, 2), M);
g(1:it) = 0;
end
