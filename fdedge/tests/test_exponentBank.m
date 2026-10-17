classdef test_exponentBank < matlab.unittest.TestCase
    % M4: exponentBank - zero-of-line recovery, zero-slope center member,
    % out-of-bank clamp.
    methods (Test)
        function zeroOfLine(t)
            % c1 = 1.01*(beta - a) with beta = 0.63 -> bhat = 0.63 (the
            % line's zero is exact; measured |bhat - 0.63| ~ 1e-16).
            % iaStar = the bank member closest to 0.63 = 0.65 (1-based
            % index 3).
            alphas = (0.55:0.05:0.85)';
            c1 = 1.01*(0.63 - alphas);
            [bhat, iaStar] = exponentBank(c1, alphas);
            t.verifyLessThan(abs(bhat - 0.63), 0.01);
            t.verifyEqual(iaStar, 3);
            t.verifyEqual(alphas(iaStar), 0.65);
        end

        function zeroSlopeCenter(t)
            % Zero-slope c1 (constant) -> p(1) = 0 -> bhat = the bank's
            % center member alphas(round(numel(alphas)/2)). For the
            % 7-member PoC bank that is alphas(4) = 0.70 (1-based
            % indexing; round(3.5) = 4) - NOT the arithmetic middle 0.75.
            alphas = (0.55:0.05:0.85)';
            c1 = 0.5*ones(size(alphas));
            [bhat, iaStar] = exponentBank(c1, alphas);
            t.verifyEqual(bhat, alphas(round(numel(alphas)/2)));
            t.verifyEqual(bhat, 0.70);
            t.verifyEqual(iaStar, round(numel(alphas)/2));
        end

        function outOfBankClamp(t)
            % Zero at beta = 1.2 (outside [0.55, 0.85]) -> clamped to the
            % top bank member 0.85.
            alphas = (0.55:0.05:0.85)';
            c1 = 1.01*(1.2 - alphas);
            [bhat, iaStar] = exponentBank(c1, alphas);
            t.verifyEqual(bhat, 0.85);
            t.verifyEqual(iaStar, numel(alphas));
        end
    end
end
