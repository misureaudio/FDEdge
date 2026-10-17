classdef test_slopePeak < matlab.unittest.TestCase
    % M3: slopePeak - localization on the PoC's smoothed pipeline, P3
    % window rejection, documented unsmoothed-ringing behavior.
    methods (Test)
        function localizesSmoothedEdge(t)
            % The PoC's actual pipeline: one-sided power-law edge,
            % Gaussian smoothing (s_grid = 8), then GL response at
            % a = 0.7 (step response since a = beta). Measured 3.55e-4
            % (0.36 grid) at M = 1024, ell = 6, window [0.40, 0.44].
            % NOTE: the plan's spec omitted the smoothing; unsmoothed the
            % GL response rings at the corner and the peak is off by
            % 4.6 grid (see unsmoothedRingingDocumented).
            M = 1024;
            x = (0:M)'/M;
            g = max(x - 0.42, 0).^0.7;
            gs = smoothGauss(g, 8);
            R = caputoGL(x, gs, 0.7);
            xhat = slopePeak(R, x, 6, 0.40, 0.44);
            t.verifyLessThan(abs(xhat - 0.42), 2/M);
        end

        function p3SpikeRejected(t)
            % P3: a huge spike at x ~= 0.99 (the spurious zero-pad
            % boundary edge, amplified) must not be selected when the
            % search is window-restricted to [0.40, 0.44]. Without the
            % window the peak jumps to ~0.98 (measured 0.9844), proving
            % the window is what rejects it.
            M = 1024;
            x = (0:M)'/M;
            g = max(x - 0.42, 0).^0.7;
            gs = smoothGauss(g, 8);
            R = caputoGL(x, gs, 0.7);
            ik0 = round(0.99*M) + 1;   % 1-based index of x ~= 0.99
            R(ik0) = R(ik0) + 1e3;
            xhat = slopePeak(R, x, 6, 0.40, 0.44);
            t.verifyLessThan(abs(xhat - 0.42), 2/M);
            xhatNw = slopePeak(R, x, 6);
            t.verifyGreaterThan(xhatNw, 0.9);
        end

        function unsmoothedRingingDocumented(t)
            % Without the PoC's Gaussian smoothing the GL response rings
            % at the corner and the peak is biased ~ell grid left of the
            % edge (measured 4.5e-3 = 4.6 grid). The PoC always smooths
            % first; this test pins the documented behavior so a future
            % "simplification" cannot silently remove the smoothing.
            M = 1024;
            x = (0:M)'/M;
            g = max(x - 0.42, 0).^0.7;
            R = caputoGL(x, g, 0.7);
            xhat = slopePeak(R, x, 6, 0.40, 0.44);
            t.verifyLessThan(abs(xhat - 0.42), 6/M);
        end
    end
end
