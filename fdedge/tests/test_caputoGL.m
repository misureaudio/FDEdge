classdef test_caputoGL < matlab.unittest.TestCase
    % M2: caputoGL - analytic reference, row/column equivalence (D13).
    methods (Test)
        function analyticT2(t)
            % D^a t^2 = 2 t^(2-a) / Gamma(3-a) exactly (Caputo, f(0)=f'(0)=0).
            % Interior point i=2048: the left-boundary layer has decayed.
            N = 2048; a = 0.7;
            tt = (0:N)'/N;
            f = tt.^2;
            D = caputoGL(tt, f, a);
            exact = 2*tt(2048)^(2-a)/gamma(3-a);
            t.verifyLessThan(abs(D(2048) - exact)/exact, 1e-3);
        end

        function rowInputEqualsColumn(t)
            % D13: row inputs give a result isequal to column inputs.
            N = 256;
            tt = (0:N)'/N;
            f = sin(2*pi*tt) + 0.1*tt;
            Dc = caputoGL(tt, f, 0.7);
            Dr = caputoGL(tt', f', 0.7);
            t.verifyTrue(isequal(Dc, Dr));
        end
    end
end
