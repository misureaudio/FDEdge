function grid = gridHybrid(cfg, xes, betas)
% Hybrid bias grids on noise-free symmetric edges, 2-pass pipeline
% (exponent at the reference position -> cusp-min position -> exponent at
% the fresh position):
%   dbeta(ie, ib) = raw 2-pass exponent bias (bcal - betas(ib))
%   dxg(ie, ia, ib) = cusp-min position bias (xhat - xes(ie))
% for the edge A| x - xes(ie) |^betas(ib), per bank order alphas(ia).
% (Move of the calibration loop in FD_EdgeDetector_PoC/
% edge_fd_hybrid_track.m; dual grids x = (0:M)'/M (Caputo/GL) and
% xr = (0:M-1)'/M (periodic, Riesz).)
% Required cfg fields: M, A, alphas, sGrid, w0, w1.
M = cfg.M;
x  = (0:M)'/M;                     % M+1-point grid (Caputo/GL)
xr = (0:M-1)'/M;                   % M-point periodic grid (Riesz)
kw = (-3*cfg.sGrid:3*cfg.sGrid)';
kg = exp(-kw.^2/(2*cfg.sGrid^2)); kg = kg/sum(kg);
dbeta = zeros(numel(xes), numel(betas));
dxg   = zeros(numel(xes), numel(cfg.alphas), numel(betas));
for ie = 1:numel(xes)
    for ib = 1:numel(betas)
        fe = cfg.A*abs(x - xes(ie)).^betas(ib);
        fs = conv(fe, kg, 'same');
        Rbank = rieszLaplacian(fs(1:M), cfg.alphas);
        [b1, ia1] = estExponent(fs, cfg.alphas, xes(ie), x, cfg.w0, cfg.w1);
        xcal = cuspMin(Rbank(:,ia1), xr, xes(ie)-0.02, xes(ie)+0.02);
        [bcal, ~] = estExponent(fs, cfg.alphas, xcal, x, cfg.w0, cfg.w1);
        dbeta(ie, ib) = bcal - betas(ib);
        for ia = 1:numel(cfg.alphas)
            dxg(ie, ia, ib) = cuspMin(Rbank(:,ia), xr, xes(ie)-0.02, xes(ie)+0.02) - xes(ie);
        end
    end
end
grid = struct('dbeta', dbeta, 'dxg', dxg);
end
