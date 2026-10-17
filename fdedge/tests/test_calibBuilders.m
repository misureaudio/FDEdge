classdef test_calibBuilders < matlab.unittest.TestCase
    % M5: the four dedicated bias-grid builders (D5) + the D7 calibrator
    % properties (A-invariance, sigma sensitivity).
    %
    % MEASURED ENVELOPE (probe_m5, deterministic):
    %   scalarHalfCross delta = 1.36943797451438698e-04
    %     -> sprintf('%.2e') = '1.37e-04' = the PoC's printed value in
    %        golden/results_poc.txt (exact string match).
    %   gridSlopePeak: [7 4], finite, deterministic; sGrid 8 vs 12 ->
    %     max abs diff = 1.512e-03 (> 1e-3: the grid IS sigma-sensitive).
    %   gridRiesz / gridHybrid on the PoC 5x7 grid: dbeta [5 7],
    %     dxg [5 7 7], finite, deterministic.
    %   A-invariance (M=512, xes=[0.3 0.4], betas=[0.6 0.7], A=1 vs 2):
    %     riesz: dbeta/dxg maxdiff = 0.0 (bit-exact); hybrid: dbeta
    %     maxdiff = 1.11e-16 (the loglogSlope 1e-6 floor), dxg = 0.0.
    methods (Test)
        function scalarHalfCrossMatchesGolden(t)
            % Finite; the %.2e print equals the PoC's printed delta in
            % golden/results_poc.txt; the full-precision value is pinned
            % to the measured reference (probe_m5).
            cfg = test_calibBuilders.baseCfg();
            delta = scalarHalfCross(cfg);
            t.verifyTrue(isfinite(delta));
            txt = fileread(golden_file('results_poc.txt'));
            tok = regexp(txt, 'estimator bias[^\n]*delta = ([-0-9.eE+]+)', 'tokens');
            t.verifyEqual(sprintf('%.2e', delta), tok{1}{1});
            t.verifyLessThan(abs(delta - 1.36943797451438698e-04), 1e-12);
        end

        function gridSlopePeakShapeDeterministic(t)
            cfg = test_calibBuilders.baseCfg();
            d = gridSlopePeak(cfg);
            t.verifyEqual(size(d), [7 4]);
            t.verifyTrue(all(isfinite(d(:))));
            t.verifyEqual(gridSlopePeak(cfg), d);
        end

        function gridRieszShapeDeterministic(t)
            cfg = test_calibBuilders.baseCfg();
            g = gridRiesz(cfg, test_calibBuilders.poCxes(), cfg.alphas);
            t.verifyEqual(size(g.dbeta), [5 7]);
            t.verifyEqual(size(g.dxg), [5 7 7]);
            t.verifyTrue(all(isfinite(g.dbeta(:))) && all(isfinite(g.dxg(:))));
            t.verifyEqual(gridRiesz(cfg, test_calibBuilders.poCxes(), cfg.alphas), g);
        end

        function gridHybridShapeDeterministic(t)
            cfg = test_calibBuilders.baseCfg();
            g = gridHybrid(cfg, test_calibBuilders.poCxes(), cfg.alphas);
            t.verifyEqual(size(g.dbeta), [5 7]);
            t.verifyEqual(size(g.dxg), [5 7 7]);
            t.verifyTrue(all(isfinite(g.dbeta(:))) && all(isfinite(g.dxg(:))));
            t.verifyEqual(gridHybrid(cfg, test_calibBuilders.poCxes(), cfg.alphas), g);
        end

        function aInvariance(t)
            % D7a: A-invariance (scale-free estimators). Tiny grid
            % (M=512, 2x2): A=1 vs A=2 agree to <= 1e-12.
            cfg = test_calibBuilders.baseCfg(); cfg.M = 512;
            xa = [0.3 0.4]'; ba = [0.6 0.7]';
            r1 = gridRiesz(cfg, xa, ba);
            cfg2 = cfg; cfg2.A = 2.0;
            r2 = gridRiesz(cfg2, xa, ba);
            t.verifyLessThan(max(abs(r1.dbeta - r2.dbeta), [], 'all'), 1e-12);
            t.verifyLessThan(max(abs(r1.dxg - r2.dxg), [], 'all'), 1e-12);
            h1 = gridHybrid(cfg, xa, ba);
            h2 = gridHybrid(cfg2, xa, ba);
            t.verifyLessThan(max(abs(h1.dbeta - h2.dbeta), [], 'all'), 1e-12);
            t.verifyLessThan(max(abs(h1.dxg - h2.dxg), [], 'all'), 1e-12);
        end

        function sigmaSensitivity(t)
            % D7b: the slope-peak grid depends on the smoothing sigma -
            % the property behind the README rule "rebuild calibration on
            % sigma change". Measured max abs diff = 1.512e-03.
            cfg8 = test_calibBuilders.baseCfg();
            cfg12 = cfg8; cfg12.sGrid = 12;
            d8 = gridSlopePeak(cfg8);
            d12 = gridSlopePeak(cfg12);
            t.verifyGreaterThan(max(abs(d8 - d12), [], 'all'), 1e-3);
        end
    end
    methods (Static)
        function cfg = baseCfg()
            % The PoC reference configuration (all five PoCs agree on
            % M, sGrid, alphas; the per-PoC extras are included so one
            % struct drives every builder).
            cfg = struct('M', 1024, 'A', 1.0, 'beta', 0.70, 'sGrid', 8, ...
                         'alphas', (0.55:0.05:0.85)', 'ell', 6, ...
                         'dpr', linspace(0.03,0.13,11), ...
                         'dprw', linspace(0.025,0.12,11), ...
                         'w0', 0.03, 'w1', 0.28);
        end
        function xes = poCxes()
            xes = (0.20:0.10:0.60)';
        end
    end
end
