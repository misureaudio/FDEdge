classdef KalmanCV
    % Constant-velocity Kalman filter: 2-state [x; v] or 4-state
    % [x; v; beta; omega].
    % predict/update are the PoC's verbatim matrix lines (edge_fd_poc.m
    % Part B 2-state; edge_fd_beta_track.m / edge_fd_riesz_track.m /
    % edge_fd_hybrid_track.m 4-state); init is r1's corrected version (D3):
    % the default P0 is n x n (v0's diag(obj.R(:))*[1;1e2] yielded a
    % vector - fatal dimension bug, review r1).
    % NOTE (deviation from the plan sketch): init/predict/update RETURN
    % obj - MATLAB classdef value semantics makes methods without a
    % return value no-ops. Call as:  kf = kf.init(x0);  kf = kf.predict();
    %  kf = kf.update(z);
    properties
        n; F; H; Q; R; state; P
    end
    methods
        function obj = KalmanCV(n, qV, qW, Rm)
            obj.n = n;
            if n == 2
                obj.F = [1 1; 0 1];
                obj.H = [1 0];
                obj.Q = qV * [1/3 1/2; 1/2 1];
            else
                obj.F = [1 1 0 0; 0 1 0 0; 0 0 1 1; 0 0 0 1];
                obj.H = [1 0 0 0; 0 0 1 0];
                obj.Q = [qV/3 qV/2 0 0; qV/2 qV 0 0; 0 0 qW/3 qW/2; 0 0 qW/2 qW];
            end
            obj.R = Rm;
        end
        function obj = init(obj, x0, P0)
            obj.state = x0;
            if nargin < 3
                % D3 fix (review r1): P0 must be n x n, not a vector.
                if obj.n == 2
                    P0 = eye(2) * obj.R(1,1);        % == PoC Part B [Rm 0; 0 Rm]
                else
                    P0 = blkdiag(obj.R(1,1)*eye(2), obj.R(2,2)*eye(2));
                end
            end
            obj.P = P0;
        end
        function obj = predict(obj)
            obj.state = obj.F * obj.state;
            obj.P = obj.F * obj.P * obj.F' + obj.Q;
        end
        function obj = update(obj, z)
            y = z - obj.H * obj.state;
            S = obj.H * obj.P * obj.H' + obj.R;
            K = obj.P * obj.H' / S;
            obj.state = obj.state + K * y;
            obj.P = (eye(obj.n) - K * obj.H) * obj.P;
        end
    end
end
