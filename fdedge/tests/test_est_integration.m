classdef test_est_integration < matlab.unittest.TestCase
    % M4: 2-pass robustness (D10) - the hybrid per-frame loop: pass 1
    % (exponent at the predicted position -> Riesz cusp position) then
    % pass 2 (exponent at the fresh position), on a noisy symmetric edge.
    methods (Test)
        function twoPassRobustness(t)
            % Noisy symmetric edge (M = 1024, xe = 0.35, beta = 0.70,
            % A = 1, sigma = 0.02, rng(0)); xguess offset by off grid
            % points, off in {+2, -2, +3}. Plan tolerance: |xraw - xe|
            % < 5/M, |bhat2 - beta| < 0.05, all finite.
            %
            % MEASURED ENVELOPE (probe_est_integration, same seed):
            %   off=+2: xraw err +1.57 grid, b1=0.6669, bhat2=0.6669
            %   off=-2: xraw err +1.57 grid, b1=0.6736, bhat2=0.6669
            %   off=+3: xraw err +1.57 grid, b1=0.6652, bhat2=0.6669
            % All three offsets converge to the SAME xraw (the loop's
            % attractor): the position is robust (1.57 grid << 5). The
            % pass-1 exponent b1 depends on the xguess offset (window
            % placement bias), but pass 2 at the fresh position is stable
            % at 0.6669: a SYSTEMATIC -0.033 bias vs beta = 0.70. That
            % bias is exactly what the PoC's calibration grid (dbeta)
            % measures and subtracts - the raw 2-pass estimate is biased
            % by design; the tolerance 0.05 covers it (margin 1.5x).
            rng(0);
            M = 1024;
            x = (0:M)'/M;
            xe = 0.35; beta = 0.70; A = 1.0; sig = 0.02;
            w0 = 0.03; w1 = 0.28;
            alphas = (0.55:0.05:0.85)';
            sgrid = 8;
            kw = (-3*sgrid:3*sgrid)';
            kg = exp(-kw.^2/(2*sgrid^2)); kg = kg/sum(kg);
            fe = A*abs(x - xe).^beta;
            fs = conv(fe + sig*randn(M+1, 1), kg, 'same');
            Rbank = rieszLaplacian(fs(1:M), alphas);
            offs = [2, -2, 3];
            for q = 1:numel(offs)
                xguess = xe + offs(q)/M;
                % pass 1: exponent at the prediction
                g1 = memorySubtract(fs, xguess, M);
                ka = min(M+1, max(2, round((xguess + w0)*M) + 1));
                kb = min(M+1, max(ka+2, round((xguess + w1)*M) + 1));
                dvec = x(ka:kb) - xguess;
                c1 = zeros(size(alphas));
                for ia = 1:numel(alphas)
                    D = caputoGL(x, g1, alphas(ia));
                    c1(ia) = loglogSlope(D(ka:kb), dvec);
                end
                [b1, ia1] = exponentBank(c1, alphas);
                % position: Riesz cusp minimum at the pass-1 order
                xraw = cuspMin(Rbank(:, ia1), x, xguess - 0.02, xguess + 0.02);
                % pass 2: exponent at the fresh position
                g2 = memorySubtract(fs, xraw, M);
                ka2 = min(M+1, max(2, round((xraw + w0)*M) + 1));
                kb2 = min(M+1, max(ka2+2, round((xraw + w1)*M) + 1));
                dvec2 = x(ka2:kb2) - xraw;
                c12 = zeros(size(alphas));
                for ia = 1:numel(alphas)
                    D = caputoGL(x, g2, alphas(ia));
                    c12(ia) = loglogSlope(D(ka2:kb2), dvec2);
                end
                [bhat2, ~] = exponentBank(c12, alphas);
                t.verifyTrue(isfinite(xraw) && isfinite(bhat2) && isfinite(b1), ...
                    sprintf('non-finite estimate at off = %d', offs(q)));
                t.verifyLessThan(abs(xraw - xe), 5/M);
                t.verifyLessThan(abs(bhat2 - beta), 0.05);
            end
        end
    end
end
