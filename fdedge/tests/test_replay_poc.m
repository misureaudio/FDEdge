classdef test_replay_poc < matlab.unittest.TestCase
    % M8: replay_poc reproduces golden/results_poc.txt through the library,
    % gated by the two-tier gate (stage 1 numeric, stage 2 byte-identical).
    % run_tests cd's to the project root before running the suite.
    methods (Test)
        function replayPoc(t)
            root = pwd;
            [pass, detail] = gate_replay(fullfile(root, 'examples', 'replay_poc.m'), ...
                'results_poc.txt', fullfile(root, 'golden', 'results_poc.txt'), ...
                fullfile(root, 'tests', 'replay_exceptions.md'));
            t.verifyTrue(pass, detail);
        end
    end
end
