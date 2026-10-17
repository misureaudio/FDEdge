classdef test_Tracker < matlab.unittest.TestCase
    % M7: Tracker orchestrator - per-profile noise-free frame at
    % (xe = 0.35, beta = 0.70), 5-frame smoke (KF engaged at frame 2),
    % gridLen validation (D4), and the Part A one-shot analyzeMatchedOrder.
    % Pre-verified in MATLAB (probe, 2026-10-16): noise-free errors
    % poc x 6.95e-06; beta x 8.38e-06 / b 1.41e-04; riesz x 2.46e-06 /
    % b 6.95e-03; hybrid x 1.31e-05 / b 5.44e-04; smoke final vhat = 0.01000
    % (all profiles); matchedOrder iBest -> alpha = 0.70, platMean err 1.9e-05.
    methods (Test)
        function pocProfile(t)
            cfg = config('Profile','poc');
            calib = scalarHalfCross(cfg);
            M = cfg.M; x = (0:M)'/M;
            % Noise-free frame at (0.35, 0.70) - within 2 grid units.
            trk = Tracker(cfg, calib);
            f = cfg.A*max(x - 0.35, 0).^0.70;
            [st, trk] = trk.update(f, NaN);
            t.verifyLessThan(abs(st.xhat - 0.35), 2/M);
            t.verifyTrue(isnan(st.bhat));   % position-only profile
            t.verifyTrue(isnan(st.vhat));   % KF not yet engaged
            % 5-frame smoke (KF engaged at frame 2).
            trk = Tracker(cfg, calib);
            for n = 1:5
                xe = 0.30 + 0.01*n;
                f = cfg.A*max(x - xe, 0).^0.70;
                [st, trk] = trk.update(f, NaN);
                t.verifyTrue(isfinite(st.xhat));
                if n == 1
                    t.verifyTrue(isnan(st.vhat));
                else
                    t.verifyTrue(isfinite(st.vhat));
                end
            end
            % D4 grid contract.
            t.verifyError(@() trk.update(f(1:500), NaN), 'fdedge:gridLen');
        end

        function betaProfile(t)
            cfg = config('Profile','beta');
            calib = gridSlopePeak(cfg);
            M = cfg.M; x = (0:M)'/M;
            trk = Tracker(cfg, calib);
            f = cfg.A*max(x - 0.35, 0).^0.70;
            [st, trk] = trk.update(f, 0.35);
            t.verifyLessThan(abs(st.xhat - 0.35), 5e-3);
            t.verifyLessThan(abs(st.bhat - 0.70), 0.03);
            t.verifyTrue(isnan(st.vhat));
            t.verifyTrue(isnan(st.ohat));
            trk = Tracker(cfg, calib);
            xguess = 0.30; rawPrev = NaN;
            for n = 1:5
                xe = 0.30 + 0.01*n;
                f = cfg.A*max(x - xe, 0).^0.70;
                [st, trk] = trk.update(f, xguess);
                t.verifyTrue(isfinite(st.xhat));
                t.verifyTrue(isfinite(st.bhat));
                if n == 1
                    xguess = st.xraw;
                else
                    xguess = st.xraw + (st.xraw - rawPrev);
                    t.verifyTrue(isfinite(st.vhat));
                    t.verifyTrue(isfinite(st.ohat));
                end
                rawPrev = st.xraw;
            end
            t.verifyError(@() trk.update(f(1:500), 0.35), 'fdedge:gridLen');
        end

        function rieszProfile(t)
            cfg = config('Profile','riesz');
            xes = (0.20:0.10:0.60)';
            calib = gridRiesz(cfg, xes, cfg.alphas);
            M = cfg.M; xr = (0:M-1)'/M;
            trk = Tracker(cfg, calib);
            f = cfg.A*abs(xr - 0.35).^0.70;
            [st, trk] = trk.update(f, 0.35);
            t.verifyLessThan(abs(st.xhat - 0.35), 5e-3);
            t.verifyLessThan(abs(st.bhat - 0.70), 0.03);
            t.verifyTrue(isnan(st.vhat));
            t.verifyTrue(isnan(st.ohat));
            trk = Tracker(cfg, calib);
            xguess = 0.30; rawPrev = NaN;
            for n = 1:5
                xe = 0.30 + 0.01*n;
                f = cfg.A*abs(xr - xe).^0.70;
                [st, trk] = trk.update(f, xguess);
                t.verifyTrue(isfinite(st.xhat));
                t.verifyTrue(isfinite(st.bhat));
                if n == 1
                    xguess = st.xraw;
                else
                    xguess = st.xraw + (st.xraw - rawPrev);
                    t.verifyTrue(isfinite(st.vhat));
                    t.verifyTrue(isfinite(st.ohat));
                end
                rawPrev = st.xraw;
            end
            t.verifyError(@() trk.update(f(1:500), 0.35), 'fdedge:gridLen');
        end

        function hybridProfile(t)
            cfg = config('Profile','hybrid');
            xes = (0.20:0.10:0.60)';
            calib = gridHybrid(cfg, xes, cfg.alphas);
            M = cfg.M; x = (0:M)'/M;
            trk = Tracker(cfg, calib);
            f = cfg.A*abs(x - 0.35).^0.70;
            [st, trk] = trk.update(f, 0.35);
            t.verifyLessThan(abs(st.xhat - 0.35), 5e-3);
            t.verifyLessThan(abs(st.bhat - 0.70), 0.03);
            t.verifyTrue(isnan(st.vhat));
            t.verifyTrue(isnan(st.ohat));
            trk = Tracker(cfg, calib);
            xguess = 0.30; rawPrev = NaN;
            for n = 1:5
                xe = 0.30 + 0.01*n;
                f = cfg.A*abs(x - xe).^0.70;
                [st, trk] = trk.update(f, xguess);
                t.verifyTrue(isfinite(st.xhat));
                t.verifyTrue(isfinite(st.bhat));
                if n == 1
                    xguess = st.xraw;
                else
                    xguess = st.xraw + (st.xraw - rawPrev);
                    t.verifyTrue(isfinite(st.vhat));
                    t.verifyTrue(isfinite(st.ohat));
                end
                rawPrev = st.xraw;
            end
            t.verifyError(@() trk.update(f(1:500), 0.35), 'fdedge:gridLen');
        end

        function matchedOrderOneShot(t)
            % Noise-free t^0.7 on the 2049-point non-normalized grid
            % t = (0:2048)'.
            cfg = config('Profile','poc');
            N = cfg.N;
            tg = (0:N)';
            f = tg.^0.7;
            res = analyzeMatchedOrder(cfg, tg, f);
            t.verifyEqual(cfg.alphas(res.iBest), 0.70);
            t.verifyLessThan(abs(res.platMean - gamma(1.7)), 1e-3);
            % ratio(iBest) is the one closest to 1 (argmin identity).
            t.verifyEqual(min(abs(res.ratio - 1)), abs(res.ratio(res.iBest) - 1));
            t.verifyEqual(res.exactPlat, cfg.A*gamma(cfg.beta + 1));
        end
    end
end
