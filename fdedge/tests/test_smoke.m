classdef test_smoke < matlab.unittest.TestCase
    methods (Test)
        function smoke(t)
            t.verifyEqual(1 + 1, 2, 'package path must load');
        end
    end
end
