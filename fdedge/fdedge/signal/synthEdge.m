function f = synthEdge(x, model, A, xe, beta, sigma)
% Synthetic edge signal: the PoC edge shapes + optional white noise.
%   model = 'oneSided'  : f = A*max(x - xe, 0).^beta
%       (poc Part A with xe = 0 on t >= 0; poc Part B; beta tracker)
%   model = 'symmetric' : f = A*abs(x - xe).^beta
%       (riesz tracker; hybrid tracker)
% REPLAY-CRITICAL RNG ORDER: when sigma > 0, EXACTLY ONE
% randn(numel(x),1) draw is added:  f = f + sigma*randn(numel(x),1).
% The replay drivers depend on this draw order matching the PoC's
% (sigB*randn(M+1,1) / sigB*randn(M,1) per frame, before smoothing).
x = x(:);
switch lower(model)
    case 'onesided'
        f = A*max(x - xe, 0).^beta;
    case 'symmetric'
        f = A*abs(x - xe).^beta;
    otherwise
        error('fdedge:synthEdge', ...
            "unknown model '%s' (expected 'oneSided' or 'symmetric')", model);
end
if sigma > 0
    f = f + sigma*randn(numel(x), 1);
end
end
