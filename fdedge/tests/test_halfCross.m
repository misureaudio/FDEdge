classdef test_halfCross < matlab.unittest.TestCase
    % M3: halfCross - step-edge localization, shift invariance (on-grid
    % exact, off-grid sub-grid envelope).
    methods (Test)
        function localizesStepEdge(t)
            % Noise-free one-sided edge at 0.55, L1 response at a = 0.7
            % (step response since a = beta), level = 0.5*quantile(R,0.9)
            % exactly as in the PoC's Part B. Measured 2.70e-4 (0.28
            % grid) at M = 1024.
            M = 1024;
            x = (0:M)'/M;
            g = max(x - 0.55, 0).^0.7;
            R = caputoL1(x, g, 0.7);
            level = 0.5*quantile(R, 0.9);
            xhat = halfCross(R, x, level);
            t.verifyLessThan(abs(xhat - 0.55), 2/M);
        end

        function shiftInvarianceOnGrid(t)
            % The property the pocB single-reference calibration relies
            % on: the bias (halfCross - xe) is constant. For edges on
            % grid points (0.25, 0.75) the L1 response is an EXACT
            % integer translation (512 grid), so the biases are equal to
            % machine precision (measured 5.6e-17).
            M = 1024;
            x = (0:M)'/M;
            xes = [0.25, 0.75];
            b = zeros(1, 2);
            for q = 1:2
                g = max(x - xes(q), 0).^0.7;
                R = caputoL1(x, g, 0.7);
                b(q) = halfCross(R, x, 0.5*quantile(R, 0.9)) - xes(q);
            end
            t.verifyLessThan(abs(b(1) - b(2)), 1e-12);
        end

        function shiftInvarianceOffGrid(t)
            % Off-grid edges (0.30, 0.70) sit at different fractional
            % cell positions (0.2 vs 0.8 of a grid cell), so the discrete
            % responses are NOT exact translations: the bias varies
            % sub-grid (measured 7.0e-5 = 0.07 grid). The plan's 1e-12
            % target is unattainable off-grid; the PoC's single-reference
            % calibration absorbs this as per-frame residual.
            M = 1024;
            x = (0:M)'/M;
            xes = [0.30, 0.70];
            b = zeros(1, 2);
            for q = 1:2
                g = max(x - xes(q), 0).^0.7;
                R = caputoL1(x, g, 0.7);
                b(q) = halfCross(R, x, 0.5*quantile(R, 0.9)) - xes(q);
            end
            t.verifyLessThan(abs(b(1) - b(2)), 1/M);
        end
    end
end
