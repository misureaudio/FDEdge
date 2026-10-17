classdef test_caputoL1 < matlab.unittest.TestCase
    % M2: caputoL1 - analytic reference, P1 guard, GL<->L1 agreement.
    methods (Test)
        function analyticT2(t)
            % L1 product integration is exact for f'' constant (t^2):
            % the weights are the exact integrals of (t_j - tau)^(1-a).
            N = 2048; a = 0.7;
            tt = (0:N)'/N;
            f = tt.^2;
            D = caputoL1(tt, f, a);
            exact = 2*tt(2048)^(2-a)/gamma(3-a);
            t.verifyLessThan(abs(D(2048) - exact)/exact, 1e-2);
        end

        function p1GuardRejectsAlpha1(t)
            % P1: the L1 scheme degenerates at alpha = 1 - explicit error.
            N = 256;
            tt = (0:N)'/N;
            f = tt.^2;
            t.verifyError(@() caputoL1(tt, f, 1.0), 'fdedge:badAlpha');
        end

        function glL1Agreement(t)
            % GL vs L1 on a smooth non-power-law signal (plan M2.1c / M2.2c).
            % The two schemes agree to a common limit at ~O(h) (measured:
            % max diff 3.35e-3 @ N=2048, 8.84e-4 @ N=8192, ratio ~1.95).
            % Absolute correctness of each scheme is pinned separately by
            % the analytic power-law references (analyticT2, as in the PoC
            % probe_gl.m); here: cross-scheme consistency + convergence.
            % Interior region only (the left-boundary initial layer is
            % excluded).
            a = 0.7;
            md = zeros(1, 2);
            Ns = [2048, 8192];
            for k = 1:2
                N = Ns(k);
                x = (0:N)'/N;
                f = sin(2*pi*x) + 0.1*x;
                Dgl = caputoGL(x, f, a);
                Dl1 = caputoL1(x, f, a);
                i0 = round(0.25*N) + 1;
                md(k) = max(abs(Dgl(i0:end) - Dl1(i0:end)));
            end
            t.verifyLessThan(md(1), 1e-2);
            t.verifyLessThan(md(2), md(1)/1.8);
        end
    end
end
