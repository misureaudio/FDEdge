classdef test_flatnessRatio < matlab.unittest.TestCase
    % M4: flatnessRatio - matched-order flatness, plateau vs Gamma(1.7).
    methods (Test)
        function matchedOrderFlatness(t)
            % The plan's spec: noise-free t^0.7 (NO smoothing), L1 bank
            % (0.40:0.05:0.95)', N = 2048, w = 128 - the PoC Part A
            % setup (t = (0:N)' in grid-index units). Measured:
            % ratio(0.70) = 0.9999 -> argmin|ratio-1| = 0.70;
            % platMean = 0.908658 vs Gamma(1.7) = 0.908639 (diff 1.9e-5
            % < 1e-3); the c1 column estimates beta - alpha to ~1e-3.
            % (With the PoC's Gaussian smoothing the plateau mean is
            % pulled down 5.8e-3 by the kernel near t = 0 - outside the
            % plan's 1e-3; the spec's "noise-free" reading wins.)
            N = 2048;
            tvec = (0:N)';
            f = tvec.^0.7;
            alphas = (0.40:0.05:0.95)';
            R = zeros(N+1, numel(alphas));
            for j = 1:numel(alphas)
                R(:, j) = caputoL1(tvec, f, alphas(j));
            end
            [~, iBeta] = min(abs(alphas - 0.70));
            [ratio, c1, platMean, iBest] = flatnessRatio(R, alphas, N, 128, iBeta);
            t.verifyEqual(alphas(iBest), 0.70);
            t.verifyLessThan(abs(platMean - gamma(1.7)), 1e-3);
            % bonus: the c1 column estimates beta - alpha (measured ~1e-3)
            t.verifyLessThan(max(abs(c1 - (0.70 - alphas))), 2e-3);
        end
    end
end
