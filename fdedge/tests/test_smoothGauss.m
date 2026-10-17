classdef test_smoothGauss < matlab.unittest.TestCase
    % M2: smoothGauss - interior constant preservation, delta mass.
    % Docstring states P3 (zero-pad boundary artifact).
    methods (Test)
        function constantPreservedInterior(t)
            % The kernel is normalized (sum(kg) = 1): a constant is
            % preserved to roundoff away from the zero-padded ends.
            N = 256;
            f = ones(N, 1);
            fs = smoothGauss(f, 8);
            interior = (3*8+2):(N-3*8-1);
            t.verifyLessThan(max(abs(fs(interior) - 1)), 1e-12);
        end

        function deltaMassInterior(t)
            % A delta's mass spreads over the 6-sigma kernel and is
            % preserved in the interior (no boundary involvement).
            N = 256;
            f = zeros(N, 1);
            f(128) = 1;
            fs = smoothGauss(f, 8);
            win = (128-3*8-2):(128+3*8+2);
            t.verifyLessThan(abs(sum(fs(win)) - 1), 1e-9);
        end
    end
end
