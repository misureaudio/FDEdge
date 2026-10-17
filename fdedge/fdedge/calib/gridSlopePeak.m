function dgrid = gridSlopePeak(cfg)
% Slope-peak bias grid: dgrid(ia, ib) = slope-peak position error of the
% GL derivative of order alphas(ia) on a noise-free one-sided edge
% (max(x - xref, 0))^beta_ref(ib) at fixed xref = 0.45, with the search
% restricted to +/-0.05 around xref (a global search would lock onto the
% spurious boundary edge).
% (Move of the calibration block in FD_EdgeDetector_PoC/
% edge_fd_beta_track.m; A = 1 as in the PoC - the grid is A-invariant by
% linearity.)
% Required cfg fields: M, alphas, sGrid, ell.
M = cfg.M;
x = (0:M)'/M;
kw = (-3*cfg.sGrid:3*cfg.sGrid)';
kg = exp(-kw.^2/(2*cfg.sGrid^2)); kg = kg/sum(kg);
beta_ref = [0.55 0.65 0.75 0.85]';
xref = 0.45;
dgrid = zeros(numel(cfg.alphas), numel(beta_ref));
for ib = 1:numel(beta_ref)
    for ia = 1:numel(cfg.alphas)
        g  = (max(x - xref, 0)).^beta_ref(ib);
        gs = conv(g, kg, 'same');
        R  = caputoGL(x, gs, cfg.alphas(ia));
        dgrid(ia, ib) = slopePeak(R, x, cfg.ell, xref - 0.05, xref + 0.05) - xref;
    end
end
end
