classdef test_replay_beta < matlab.unittest.TestCase
    % M9: replay_beta reproduces golden/results_beta.txt through the
    % library, gated by the two-tier gate (stage 1 numeric, stage 2
    % byte-identical). run_tests cd's to the project root before running
    % the suite.
    methods (Test)
        function replayBeta(t)
            root = pwd;
            [pass, detail] = gate_replay(fullfile(root, 'examples', 'replay_beta.m'), ...
                'results_beta.txt', fullfile(root, 'golden', 'results_beta.txt'), ...
                fullfile(root, 'tests', 'replay_exceptions.md'));
            t.verifyTrue(pass, detail);
        end
    end
end
