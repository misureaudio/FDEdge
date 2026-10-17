function grid = gridRiesz(cfg, xes, betas)
% Riesz bias grids on noise-free symmetric edges:
%   dbeta(ie, ib) = raw Riesz exponent bias (braw - betas(ib))
%   dxg(ie, ia, ib) = cusp-min position bias (xhat - xes(ie))
% for the edge A| x - xes(ie) |^betas(ib), per bank order alphas(ia).
% (Move of the calibration loop in FD_EdgeDetector_PoC/
% edge_fd_riesz_track.m; the periodic M-point grid x = (0:M-1)'/M is the
% Riesz PoC convention.)
% Required cfg fields: M, A, alphas, sGrid, dpr, dprw.
M = cfg.M;
x = (0:M-1)'/M;                     % periodic grid (Riesz PoC convention)
kw = (-3*cfg.sGrid:3*cfg.sGrid)';
kg = exp(-kw.^2/(2*cfg.sGrid^2)); kg = kg/sum(kg);
dbeta = zeros(numel(xes), numel(betas));
dxg   = zeros(numel(xes), numel(cfg.alphas), numel(betas));
for ie = 1:numel(xes)
    for ib = 1:numel(betas)
        fs = conv(cfg.A*abs(x - xes(ie)).^betas(ib), kg, 'same');
        Rbank = rieszLaplacian(fs, cfg.alphas);
        [braw, ~, ~] = estimateFrame(Rbank, cfg.alphas, xes(ie), x, cfg.dpr, cfg.dprw);
        dbeta(ie, ib) = braw - betas(ib);
        for ia = 1:numel(cfg.alphas)
            dxg(ie, ia, ib) = cuspMin(Rbank(:,ia), x, xes(ie)-0.02, xes(ie)+0.02) - xes(ie);
        end
    end
end
grid = struct('dbeta', dbeta, 'dxg', dxg);
end

function [bhat, xhat, iaStar] = estimateFrame(Rbank, alphas, xref, xg, dpr, dprw)
% Per-frame Riesz estimator: gfit exponents (model R = C0 + C*d^g) +
% cusp-min position. (Move of the local function in
% FD_EdgeDetector_PoC/edge_fd_riesz_track.m, with the package
% rieszGfit / cuspMin / exponentBank inlined.)
M = numel(xg);
i0  = round(xref*M);
ip  = round(dpr *M);
ipw = round(dprw*M);
% right-side probes, modular indices (the periodic response depends only
% on the circular distance, so wrapped probes are at the stated distance).
iR  = mod(i0 + ip  - 1, M) + 1;
iRw = mod(i0 + ipw - 1, M) + 1;
ghat = zeros(size(alphas));
for ia = 1:numel(alphas)
    [gh, ok] = rieszGfit(Rbank(iR, ia), dpr);
    if ~ok
        [gh, ok] = rieszGfit(Rbank(iRw, ia), dprw);
    end
    if ~ok
        gh = 0;
    end
    ghat(ia) = gh;
end
% beta_hat = zero of the (a, ghat(a)) line (ghat(a) ~ beta - a + bias(a))
[bhat, iaStar] = exponentBank(ghat, alphas);
xhat = cuspMin(Rbank(:,iaStar), xg, xref - 0.02, xref + 0.02);
end
