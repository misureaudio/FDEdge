function nums = gate_parse(text)
% Extract every printed number from a results text with a tolerance of
% 2x the PLACE VALUE of its last printed digit (printed integers: exact).
%   '0.2008'    -> v=0.2008,   tol=2e-4   (last digit at 1e-4 place)
%   '1.932e-03' -> v=1.932e-3, tol=2e-6   (last digit at 1e-6 place)
%   '1024'      -> v=1024,     tol=0      (printed integer: exact)
%   '-0.1'      -> v=-0.1,     tol=2e-1
% This is the tolerance rule of the two-tier replay gate (stage 1):
% portable across platforms (last-ulp flips that move a printed digit
% are tolerated), strict enough to catch any real regression.
% NOTE: plain capturing groups + 'tokens' output are used because the
% named-group syntax (?P<name>...) silently fails to match in this
% MATLAB R2024a environment (verified by probe).
C = regexp(text, '[-+]?(?:\d+\.\d+|\d+\.|\d+)(?:e[+-]\d+)?', 'match');
nums = struct('v', {}, 'tol', {});
for i = 1:numel(C)
    s = C{i};
    sign = 1;
    if s(1) == '-' || s(1) == '+'
        if s(1) == '-'; sign = -1; end
        s = s(2:end);
    end
    tok = regexp(s, '^(\d+\.\d+|\d+\.|\d+)(e[+-]\d+)?$', 'tokens');
    if isempty(tok)
        continue;                                   % defensive: cannot happen for outer matches
    end
    parts = tok{1};                                 % {mantissa, exponent}
    mant = parts{1};
    e = 0;
    hasExp = ~isempty(parts{2});
    if hasExp
        e = str2double(parts{2}(2:end));
    end
    if contains(mant, '.')
        k = numel(mant) - 1;                        % printed decimal digits
        tol = 2 * 10^(e - k);
    elseif hasExp
        tol = 2 * 10^(e - (numel(mant) - 1));       % integer mantissa + exponent
    else
        tol = 0;                                    % printed integer: exact
    end
    nums(end+1) = struct('v', sign * str2double(s), 'tol', tol);
end
end
