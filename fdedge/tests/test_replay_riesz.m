classdef test_replay_riesz < matlab.unittest.TestCase
    % M10: replay_riesz reproduces golden/results_riesz.txt through the
    % library, gated by the two-tier gate (stage 1 numeric, stage 2
    % byte-identical). run_tests cd's to the project root before running
    % the suite.
    methods (Test)
        function replayRiesz(t)
            root = pwd;
            [pass, detail] = gate_replay(fullfile(root, 'examples', 'replay_riesz.m'), ...
                'results_riesz.txt', fullfile(root, 'golden', 'results_riesz.txt'), ...
                fullfile(root, 'tests', 'replay_exceptions.md'));
            t.verifyTrue(pass, detail);
        end
    end
end
