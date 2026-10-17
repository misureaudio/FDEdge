classdef test_rieszLaplacian < matlab.unittest.TestCase
    % M2: rieszLaplacian - analytic Gaussian reference, P4 shape checks,
    % M-limit guard.
    methods (Test)
        function gaussianReference(t)
            % (-Delta)^{a/2} of the unit-mass Gaussian C*exp(-x^2/4),
            % C = (4*pi)^(-1/2), at x = 0. Continuum value (exact): the
            % Gaussian's FT is exp(-k^2), so the value is
            % (1/2pi)*int |k|^0.7 exp(-k^2) dk = gamma(0.85)/(2*pi)
            % = 0.1770572858.
            % The windowed-FFT operator computes the PERIOD-L operator
            % (L = M*dx). Because the |k|^a multiplier has a cusp at k=0
            % (C^0, not C^1), the trapezoidal rule over the k grid carries
            % algebraic aliasing ~L^{-1-a} (measured: 9.60e-3 @ L=16,
            % error ratio 3.250 per 2x L ~= 2^1.7). Convergence to the
            % continuum value requires L -> inf; refining dx alone at fixed
            % L only re-samples the same period-L operator (M=1024 and
            % M=4096 at L=16 give identical values to 1e-10).
            a = 0.7;
            dx = 1/64;
            Ls = [64, 128];
            errs = zeros(1, 2);
            for q = 1:2
                L = Ls(q);
                M = L/dx;
                x = (-L/2 + (0:M-1)'*dx);
                f = (4*pi)^(-0.5)*exp(-x.^2/4);
                R = rieszLaplacian(f, a, dx);
                errs(q) = R(M/2+1) - gamma(0.85)/(2*pi);
            end
            % L=128 (M=8192, sandbox-safe limit): measured 2.77e-4.
            t.verifyLessThan(abs(errs(2)), 5e-4);
            % convergence rate: measured 3.25 ~= 2^1.7 per doubling of L
            t.verifyGreaterThan(abs(errs(1))/abs(errs(2)), 2.5);
        end

        function p4ShapeChecks(t)
            % P4: scalar order -> vector (a row k would make M x M -> OOM);
            % bank -> M x 7.
            M = 1024;
            f = sin(2*pi*(0:M-1)'/M);
            R = rieszLaplacian(f, 0.7);
            t.verifyEqual(size(R), [M 1]);
            Rb = rieszLaplacian(f, (0.55:0.05:0.85)');
            t.verifyEqual(size(Rb), [M 7]);
        end

        function bigMGuard(t)
            % M > 8192 is refused before any FFT (sandbox memory limit).
            f = sin(2*pi*(0:16383)'/16384);
            t.verifyError(@() rieszLaplacian(f, 0.7), 'fdedge:bigM');
        end
    end
end
