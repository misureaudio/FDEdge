function [ghat, ok] = rieszGfit(y, d)
% Fit R(d) = C*d^g + C0 (3-parameter) by bounded grid search over g
% (for fixed g the (C, C0) fit is linear) + parabolic refinement.
% VALIDATED (probe_est_variants): over 11 probes on [0.03, 0.13] the
% noise-free g error is <= 0.04 (structured in xref, beta -> absorbed by
% the calibration grid). A quadratic background (5-param) destroys g
% identifiability (drifts to the grid edge) - do NOT add b1, b2. A wider
% range (to 0.20) or a nearer start (0.02) biases g by -0.15..-0.25.
% Returns ok=false when the residual surface is flat (no resolvable
% cusp) or the inputs are degenerate.
% P5 (documented limitation): under noise near a domain boundary the
% 3-parameter fit runs AWAY to the grid edge (the biased minimum settles
% at the range boundary). Do NOT narrow the range, do NOT add quadratic
% background terms, do NOT symmetrize - the PoC's calibration grid
% absorbs the structured bias.
% Moved verbatim from FD_EdgeDetector_PoC/edge_fd_riesz_track.m: gfit
% (output renamed g -> ghat).
y = y(:); d = d(:);   % column orientation (row inputs from the caller)
if ~all(isfinite(y))
    ghat = 0; ok = false;
    return;
end
% g range [-0.45, 0.45]: the true g = beta-a is in [-0.30, +0.30], but the
% curved periodic-extension background biases the (C,g,C0) fit toward the
% range edge; a wider range lets the fit settle at its (biased) minimum
% instead of being clipped, which gives a better line fit. (A +/-0.35
% range clips the biased minimum and degrades the exponent RMSE.)
gt = (-0.45:0.0225:0.45)';
Rres = zeros(size(gt));
for j = 1:numel(gt)
    if abs(gt(j)) < 1e-12
        A = [ones(numel(d),1)];   % g = 0: the d^g column degenerates to a constant
    else
        A = [d.^gt(j) ones(numel(d),1)];
    end
    p = A \ y;
    Rres(j) = sum((y - A*p).^2);
end
if ~all(isfinite(Rres))
    ghat = 0; ok = false;
    return;
end
[rmin, jmin] = min(Rres);   % FIRST output = value, second = index
% flatness test relative to the noise floor (rmin ~ (n-p)*sigma^2)
if (max(Rres) - rmin) < 2.5*rmin
    ghat = 0; ok = false;   % flat residual surface: g not resolvable
    return;
end
if jmin > 1 && jmin < numel(gt)
    d2 = Rres(jmin-1) - 2*Rres(jmin) + Rres(jmin+1);
    if abs(d2) > 1e-14
        ghat = gt(jmin) + 0.5*(Rres(jmin-1) - Rres(jmin+1))/d2 * (gt(2)-gt(1));
    else
        ghat = gt(jmin);
    end
else
    ghat = gt(jmin);
end
ghat = min(max(ghat, -0.45), 0.45);
ok = true;
end
