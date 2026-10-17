function plotTracks(frames, xe, xhatRaw, xhatKf, betat, bhatRaw, bhatKf, outPng, titleSuf)
%PLOTTRACKS Two-panel tracking figure, PoC partD/partE style.
%   Left panel:  edge position vs frame - true (black solid), raw
%   debiased (red dots), Kalman (blue dashed).
%   Right panel: edge exponent vs frame - same three traces.
%   PLOTTRACKS(frames, xe, xhatRaw, xhatKf, betat, bhatRaw, bhatKf,
%   outPng, titleSuf) exports the figure to OUTPNG at 150 dpi and
%   closes it.
%
% Visual only - never gated. The replay drivers skip PNG output; this
% serves examples/demo.m (demo.png).
fig = figure('visible','off');
subplot(1,2,1);
plot(frames, xe, 'k-', 'LineWidth', 1.6); hold on;
plot(frames, xhatRaw, 'r.', 'MarkerSize', 14);
plot(frames, xhatKf, 'b--', 'LineWidth', 1.4);
xlabel('frame'); ylabel('edge position x_e');
legend('true', 'raw (debiased)', 'Kalman', 'Location', 'southoutside');
title(['Position tracking' titleSuf]);
subplot(1,2,2);
plot(frames, betat, 'k-', 'LineWidth', 1.6); hold on;
plot(frames, bhatRaw, 'r.', 'MarkerSize', 14);
plot(frames, bhatKf, 'b--', 'LineWidth', 1.4);
xlabel('frame'); ylabel('edge exponent \beta');
legend('true', 'raw \beta_hat', 'Kalman', 'Location', 'southoutside');
title(['Exponent tracking' titleSuf]);
exportgraphics(fig, outPng, 'Resolution', 150);
close(fig);
end
