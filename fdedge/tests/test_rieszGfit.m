classdef test_rieszGfit < matlab.unittest.TestCase
    % M4: rieszGfit - 3-parameter grid-search fit, flatness rejection,
    % degenerate-input guard.
    methods (Test)
        function noiseFreeFit(t)
            % C*d^g + C0 with g = 0.12, C = 1.5, C0 = -0.3 on the PoC's
            % probe grid d = [0.03, 0.13] (11 pts). Measured ghat =
            % 0.120010 (grid + parabolic refinement), ok = true.
            d = linspace(0.03, 0.13, 11);
            y = 1.5*d.^0.12 - 0.3;
            [ghat, ok] = rieszGfit(y, d);
            t.verifyTrue(ok);
            t.verifyLessThan(abs(ghat - 0.12), 0.01);
        end

        function flatResponseRejected(t)
            % No power-law component (constant + white noise): the
            % residual surface is flat -> ok = false, ghat = 0.
            % NOTE: the PoC's flatness test (max(Rres) - rmin) < 2.5*rmin
            % only triggers when the noise floor rmin ~ (n-p)*sigma^2 is
            % nonzero; a NOISE-FREE constant gives rmin ~ 0 (exact fit at
            % every g) and slips past the test (a known degenerate case
            % of the PoC code), so this test uses noise.
            rng(42);
            d = linspace(0.03, 0.13, 11);
            y = 0.5 + 0.02*randn(11, 1);
            [ghat, ok] = rieszGfit(y, d);
            t.verifyFalse(ok);
            t.verifyEqual(ghat, 0);
        end

        function degenerateInputs(t)
            % Non-finite input -> ok = false, ghat = 0 (the PoC's guard).
            d = linspace(0.03, 0.13, 11);
            y = 1.5*d.^0.12;
            y(3) = NaN;
            [ghat, ok] = rieszGfit(y, d);
            t.verifyFalse(ok);
            t.verifyEqual(ghat, 0);
        end
    end
end
