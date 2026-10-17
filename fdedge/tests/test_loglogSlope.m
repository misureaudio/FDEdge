classdef test_loglogSlope < matlab.unittest.TestCase
    % M4: loglogSlope - power-law exactness, 1e-6 floor on zero input.
    methods (Test)
        function purePowerLaw(t)
            % C*d^g -> slope = g (polyfit of the exact log-log line).
            % Measured 4.4e-16 at g = 1.3.
            d = linspace(0.03, 0.13, 50);
            g = 1.3;
            y = 2.5*d.^g;
            c1 = loglogSlope(y, d);
            t.verifyLessThan(abs(c1 - g), 1e-10);
        end

        function allZeroInputFinite(t)
            % The 1e-6 floor keeps all-zero responses finite (slope ~ 0;
            % measured -1.6e-15).
            d = linspace(0.03, 0.13, 50);
            c1 = loglogSlope(zeros(50, 1), d);
            t.verifyTrue(isfinite(c1));
            t.verifyLessThan(abs(c1), 1e-6);
        end
    end
end
