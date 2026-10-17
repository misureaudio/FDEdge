function res = analyzeMatchedOrder(cfg, t, f)
% Part A one-shot analysis (edge_fd_poc.m Part A, lines 38-57): the L1
% response bank over the non-normalized grid t = (0:N)' plus the flatness
% metrics. Part A is NOT a tracking profile (D4) - this is the one-shot
% matched-order analysis for the 2049-point grid.
%   t: (N+1)-point grid, t = (0:N)' (non-normalized, as in the PoC)
%   f: signal on t (smoothed by the caller, if desired - the PoC smooths
%      before the bank)
% Returns struct('ratio',...,'c1',...,'iBest',...,'platMean',...,'exactPlat',...)
% (move of the PoC Part A block: fdL1 -> caputoL1, ratio/plateau block ->
% flatnessRatio; same operation order.)
N = numel(t) - 1;
R = zeros(N+1, numel(cfg.alphas));
for j = 1:numel(cfg.alphas)
    R(:,j) = caputoL1(t, f, cfg.alphas(j));
end
[~, iBeta] = min(abs(cfg.alphas - cfg.beta));
w = 128;   % PoC Part A endpoint window width (grid points)
[ratio, c1, platMean, iBest] = flatnessRatio(R, cfg.alphas, N, w, iBeta);
res = struct('ratio', ratio, 'c1', c1, 'iBest', iBest, ...
             'platMean', platMean, 'exactPlat', cfg.A*gamma(cfg.beta+1));
end
