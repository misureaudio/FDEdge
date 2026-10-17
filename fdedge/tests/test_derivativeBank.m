classdef test_derivativeBank < matlab.unittest.TestCase
    % M2: derivativeBank - dispatch facade equals the direct calls.
    methods (Test)
        function glColumnsEqualDirect(t)
            N = 256;
            tt = (0:N)'/N;
            f = sin(2*pi*tt) + 0.1*tt;
            a = (0.55:0.05:0.85)';
            R = derivativeBank(tt, f, a, 'Method', 'gl');
            t.verifyEqual(size(R), [N+1 7]);
            for ia = 1:numel(a)
                t.verifyTrue(isequal(R(:,ia), caputoGL(tt, f, a(ia))));
            end
        end

        function l1ColumnsEqualDirect(t)
            N = 256;
            tt = (0:N)'/N;
            f = sin(2*pi*tt) + 0.1*tt;
            a = (0.55:0.05:0.85)';
            R = derivativeBank(tt, f, a, 'Method', 'l1');
            for ia = 1:numel(a)
                t.verifyTrue(isequal(R(:,ia), caputoL1(tt, f, a(ia))));
            end
        end

        function rieszBankEqualsDirect(t)
            N = 256;
            tt = (0:N)'/N;
            f = sin(2*pi*tt);
            a = (0.55:0.05:0.85)';
            R = derivativeBank(tt, f, a, 'Method', 'riesz');
            t.verifyTrue(isequal(R, rieszLaplacian(f, a)));
        end
    end
end
