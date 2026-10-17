function R = derivativeBank(t, f, a, varargin)
% Facade: response bank of an FD operator over a list of orders.
%   R = derivativeBank(t, f, alphas, 'Method', 'gl'|'l1'|'riesz')
% 'gl' and 'l1' need the grid t; 'riesz' ignores t (periodic grid).
% Moved from the per-frame bank loops in the PoC scripts (same operation
% order: one operator call per order, in order - standing rule D14).
method = 'gl';
if nargin >= 4
    if mod(numel(varargin), 2) ~= 0 || ~ischar(varargin{1})
        error('fdedge:badArg', 'derivativeBank: optional arguments must be name/value pairs');
    end
    for k = 1:2:numel(varargin)
        if ~ischar(varargin{k}) || ~ischar(varargin{k+1})
            error('fdedge:badArg', 'derivativeBank: optional arguments must be name/value pairs');
        end
        if ~strcmpi(varargin{k}, 'Method')
            error('fdedge:badArg', 'derivativeBank: unknown option ''%s''', varargin{k});
        end
        method = lower(varargin{k+1});
    end
    if ~any(strcmp(method, {'gl', 'l1', 'riesz'}))
        error('fdedge:badMethod', 'derivativeBank: Method must be ''gl'', ''l1'' or ''riesz''; got ''%s''', method);
    end
end
f = f(:);
if strcmp(method, 'riesz')
    R = rieszLaplacian(f, a);
    return
end
if isscalar(a)
    if strcmp(method, 'gl')
        R = caputoGL(t, f, a);
    else
        R = caputoL1(t, f, a);
    end
    return
end
R = zeros(numel(f), numel(a));
for ia = 1:numel(a)
    if strcmp(method, 'gl')
        R(:,ia) = caputoGL(t, f, a(ia));
    else
        R(:,ia) = caputoL1(t, f, a(ia));
    end
end
end
