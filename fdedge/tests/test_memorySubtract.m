classdef test_memorySubtract < matlab.unittest.TestCase
    % M4: memorySubtract - P6 identity (GL linearity) and clean-slope
    % recovery on the PoC's smoothed pipeline.
    methods (Test)
        function p6IdentityGL(t)
            % P6: D^a f - D^a f_left == D^a g. The GL scheme is linear
            % and f = g + f_left pointwise (g is f with the left side
            % zeroed, f_left is f with the right side zeroed), so the
            % identity is exact up to the rounding of the separate
            % h^(-a)*conv sums: measured 2.8e-14, plan tolerance 1e-12.
            M = 1024;
            x = (0:M)'/M;
            xe = 0.4;
            a = 0.7;
            f = abs(x - xe).^0.7;
            g = memorySubtract(f, xe, M);
            it = min(max(round(xe*M) + 1, 2), M);
            fleft = f; fleft(it+1:end) = 0;
            diff = caputoGL(x, f, a) - caputoGL(x, fleft, a) - caputoGL(x, g, a);
            t.verifyLessThan(max(abs(diff)), 1e-12);
        end

        function cleanSlopeRecovery(t)
            % The probe_hybrid_memsub property: the log-log slope of the
            % memory-subtracted GL response over [xe+w0, xe+w1] is
            % c1 ~= 1.01*(beta - a). Measured ratio 1.0153 with the
            % PoC's Gaussian smoothing (s_grid = 8); 1.0029 unsmoothed.
            % Plan tolerance: 5%.
            M = 1024;
            x = (0:M)'/M;
            xe = 0.35; beta = 0.7; a = 0.55;
            w0 = 0.03; w1 = 0.28;
            fe = abs(x - xe).^beta;
            fs = smoothGauss(fe, 8);
            g = memorySubtract(fs, xe, M);
            ka = min(M+1, max(2, round((xe + w0)*M) + 1));
            kb = min(M+1, max(ka+2, round((xe + w1)*M) + 1));
            dvec = x(ka:kb) - xe;
            D = caputoGL(x, g, a);
            c1 = loglogSlope(D(ka:kb), dvec);
            t.verifyLessThan(abs(c1 - (beta - a))/abs(beta - a), 0.05);
        end
    end
end
