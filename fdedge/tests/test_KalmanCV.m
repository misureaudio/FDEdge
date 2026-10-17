classdef test_KalmanCV < matlab.unittest.TestCase
    % M6: KalmanCV - D3 P0 dimension fix, 2-state line tracking, 4-state
    % PoC-trajectory tracking with the PoC's explicit P0.
    % Pre-verified in MATLAB (probe, 2026-10-16): 2-state noise-free
    % v err = 4.3e-17; 4-state noise-free v/w err = 3.5e-18 / 2.3e-17;
    % 4-state noisy (rng(0), sx=1e-3, sb=1e-2) v/w err = 1.9e-05 / 1.0e-03.
    methods (Test)
        function p0Dimensions(t)
            % D3: the default P0 must be n x n (v0's diag(obj.R(:))*[1;1e2]
            % yielded a vector - fatal).
            kf2 = KalmanCV(2, 1e-7, 1e-12, [2e-7]);
            kf2 = kf2.init([0.25; 0.01]);
            t.verifyEqual(size(kf2.P), [2 2]);
            kf4 = KalmanCV(4, 1e-8, 1e-12, [2.1e-6 0; 0 1.2e-2]);
            kf4 = kf4.init([0.20; 0.01; 0.55; 0.30/39]);
            t.verifyEqual(size(kf4.P), [4 4]);
        end

        function twoStateNoiseFreeLine(t)
            % z = 0.25 + 0.01*n (n = 1..60): the CV model is exact, so the
            % velocity estimate must converge to 0.01.
            n = 1:60;
            z = 0.25 + 0.01*n;
            kf = KalmanCV(2, 1e-7, 1e-12, [1e-6]);
            kf = kf.init([z(1); z(2) - z(1)]);
            for k = 2:60
                kf = kf.predict();
                kf = kf.update(z(k));
            end
            t.verifyLessThan(abs(kf.state(2) - 0.01), 1e-4);
        end

        function fourStateNoiseFreePoCTrajectory(t)
            % The PoC trajectory (xe = 0.20 + 0.01*n', beta = 0.55 +
            % (0.30/39)*n') with the PoC's explicit P0
            % [Rx 0 0 0; 0 Rx*1e2 0 0; 0 0 Rb 0; 0 0 0 Rb*1e4]; Rx/Rb are
            % the golden raw variances (RMSE^2 from results_riesz.txt:
            % position raw RMSE 1.460e-03, exponent raw RMSE 1.098e-01).
            nfr = 40;
            xe = 0.20 + 0.01*(0:nfr-1)';
            betat = 0.55 + (0.30/(nfr-1))*(0:nfr-1)';
            Rx = (1.460e-3)^2;
            Rb = (1.098e-1)^2;
            kf = KalmanCV(4, 1e-8, 1e-12, [Rx 0; 0 Rb]);
            kf = kf.init([xe(1); xe(2)-xe(1); betat(1); betat(2)-betat(1)], ...
                [Rx 0 0 0; 0 Rx*1e2 0 0; 0 0 Rb 0; 0 0 0 Rb*1e4]);
            for k = 2:nfr
                kf = kf.predict();
                kf = kf.update([xe(k); betat(k)]);
            end
            t.verifyLessThan(abs(kf.state(2) - 0.01), 1e-4);
            t.verifyLessThan(abs(kf.state(4) - 0.30/39), 5e-4);
        end

        function fourStateNoisyPoCTrajectory(t)
            % Noisy variant: same KF, measurement noise sx = 1e-3,
            % sb = 1e-2 (rng(0)). Tolerances are the probe-measured
            % envelope (v: 1.9e-05 -> 1e-4; w: 1.0e-03 -> 2.5e-3) - the
            % plan's 5e-4 applies to the noise-free case above.
            rng(0);
            nfr = 40;
            xe = 0.20 + 0.01*(0:nfr-1)';
            betat = 0.55 + (0.30/(nfr-1))*(0:nfr-1)';
            zn = [xe + 1e-3*randn(nfr,1), betat + 1e-2*randn(nfr,1)];
            Rx = (1.460e-3)^2;
            Rb = (1.098e-1)^2;
            kf = KalmanCV(4, 1e-8, 1e-12, [Rx 0; 0 Rb]);
            kf = kf.init([zn(1,1); zn(2,1)-zn(1,1); zn(1,2); zn(2,2)-zn(1,2)], ...
                [Rx 0 0 0; 0 Rx*1e2 0 0; 0 0 Rb 0; 0 0 0 Rb*1e4]);
            for k = 2:nfr
                kf = kf.predict();
                kf = kf.update(zn(k,:));
            end
            t.verifyLessThan(abs(kf.state(2) - 0.01), 1e-4);
            t.verifyLessThan(abs(kf.state(4) - 0.30/39), 2.5e-3);
        end
    end
end
