classdef test_config < matlab.unittest.TestCase
    % M6: config - every profile's defaults reproduce the PoC constants
    % (the v0 test list), and the D9 winScale invariant
    % (w0 > 3*sGrid/M for beta/riesz/hybrid).
    methods (Test)
        function pocDefaults(t)
            cfg = config('Profile','poc');
            t.verifyEqual(cfg.M, 1024);
            t.verifyEqual(cfg.N, 2048);
            t.verifyEqual(cfg.alphas, (0.40:0.05:0.95)');
            t.verifyEqual(cfg.beta, 0.7);
            t.verifyEqual(cfg.sigA, 0.01);
            t.verifyEqual(cfg.sigB, 0.02);
            t.verifyEqual(cfg.qV, 1e-7);
            t.verifyEqual(cfg.nfr, 60);
            t.verifyEqual(cfg.xe0, 0.25);
            t.verifyEqual(cfg.v, 0.01);
            t.verifyEqual(cfg.A, 1.0);
            t.verifyEqual(cfg.sGrid, 8);
        end

        function trackingDefaults(t)
            % v0 test list (beta profile; riesz/hybrid identical here).
            cfg = config('Profile','beta');
            t.verifyEqual(cfg.M, 1024);
            t.verifyEqual(cfg.alphas, (0.55:0.05:0.85)');
            t.verifyEqual(cfg.w0, 0.03);
            t.verifyEqual(cfg.w1, 0.28);
            t.verifyEqual(cfg.sGrid, 8);
            t.verifyEqual(cfg.sigB, 0.02);
            t.verifyEqual(cfg.qV, 1e-8);
            t.verifyEqual(cfg.qW, 1e-12);
            t.verifyEqual(cfg.xLo, 0.17);
            t.verifyEqual(cfg.xHi, 0.83);
            t.verifyEqual(cfg.win, 0.02);
            t.verifyEqual(cfg.nfr, 40);
            t.verifyEqual(cfg.xe0, 0.20);
            t.verifyEqual(cfg.v, 0.01);
            t.verifyEqual(cfg.b0, 0.55);
            t.verifyEqual(cfg.w, 0.30/39);
            t.verifyEqual(cfg.ell, 6);
        end

        function rieszExtras(t)
            cfg = config('Profile','riesz');
            t.verifyEqual(cfg.dpr, linspace(0.03,0.13,11));
            t.verifyEqual(cfg.dprw, linspace(0.025,0.12,11));
        end

        function hybridDefaultPassesWinScale(t)
            % 0.03 > 3*8/1024 = 0.0234375 - no error.
            cfg = config('Profile','hybrid');
            t.verifyGreaterThan(cfg.w0, 3*cfg.sGrid/cfg.M);
        end

        function winScaleError(t)
            t.verifyError(@() config('Profile','hybrid','w0',0.01), ...
                'fdedge:winScale');
        end

        function overrideApplies(t)
            cfg = config('Profile','hybrid','A',2.0,'nfr',50);
            t.verifyEqual(cfg.A, 2.0);
            t.verifyEqual(cfg.nfr, 50);
            t.verifyEqual(cfg.M, 1024);   % untouched defaults remain
        end

        function unknownProfile(t)
            t.verifyError(@() config('Profile','wiggly'), 'fdedge:config');
        end

        function unknownField(t)
            t.verifyError(@() config('Profile','hybrid','bogus',1), ...
                'fdedge:config');
        end
    end
end
