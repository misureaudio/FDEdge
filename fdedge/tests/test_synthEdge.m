classdef test_synthEdge < matlab.unittest.TestCase
    % M6: synthEdge - the PoC edge shapes, and the replay-critical RNG
    % order (exactly one randn(numel(x),1) draw, added last).
    methods (Test)
        function noiseFreeShapes(t)
            x = (0:100)'/100;
            A = 2.0; xe = 0.3; beta = 0.7;
            t.verifyEqual(synthEdge(x, 'oneSided', A, xe, beta, 0), ...
                A*max(x - xe, 0).^beta);
            t.verifyEqual(synthEdge(x, 'symmetric', A, xe, beta, 0), ...
                A*abs(x - xe).^beta);
        end

        function oneSidedZeroLeft(t)
            % The one-sided shape is exactly zero at and left of the edge.
            x = (0:100)'/100;
            f = synthEdge(x, 'oneSided', 1.0, 0.3, 0.7, 0);
            t.verifyEqual(f(x <= 0.3), zeros(numel(x(x <= 0.3)), 1));
        end

        function singleRandnDraw(t)
            % REPLAY-CRITICAL: with a fresh seed, the noise must be
            % bit-exact equal to the next randn block (exactly one draw,
            % added last). Pre-verified: isequal = 1 for both models.
            x = (0:99)'/100;
            rng(0);
            f1 = synthEdge(x, 'oneSided', 1.0, 0.3, 0.7, 0.02);
            rng(0);
            f2 = 1.0*max(x - 0.3, 0).^0.7 + 0.02*randn(numel(x), 1);
            t.verifyEqual(f1, f2);
            rng(0);
            g1 = synthEdge(x, 'symmetric', 3.5, 0.42, 0.9, 0.05);
            rng(0);
            g2 = 3.5*abs(x - 0.42).^0.9 + 0.05*randn(numel(x), 1);
            t.verifyEqual(g1, g2);
        end

        function rngReproducibility(t)
            x = (0:99)'/100;
            rng(0);
            f1 = synthEdge(x, 'oneSided', 1.0, 0.3, 0.7, 0.02);
            rng(0);
            f2 = synthEdge(x, 'oneSided', 1.0, 0.3, 0.7, 0.02);
            t.verifyEqual(f1, f2);
        end

        function unknownModel(t)
            x = (0:10)'/10;
            t.verifyError(@() synthEdge(x, 'wiggly', 1, 0.5, 0.7, 0), ...
                'fdedge:synthEdge');
        end
    end
end
