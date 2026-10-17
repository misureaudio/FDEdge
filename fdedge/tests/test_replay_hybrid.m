classdef test_replay_hybrid < matlab.unittest.TestCase
    % M11: replay_hybrid reproduces golden/results_hybrid.txt through the
    % library, gated by the two-tier gate (stage 1 numeric, stage 2
    % byte-identical). run_tests cd's to the project root before running
    % the suite.
    methods (Test)
        function replayHybrid(t)
            root = pwd;
            [pass, detail] = gate_replay(fullfile(root, 'examples', 'replay_hybrid.m'), ...
                'results_hybrid.txt', fullfile(root, 'golden', 'results_hybrid.txt'), ...
                fullfile(root, 'tests', 'replay_exceptions.md'));
            t.verifyTrue(pass, detail);
        end
    end
end
