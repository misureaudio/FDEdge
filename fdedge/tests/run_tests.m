%% run_tests.m - headless unittest runner, exits non-zero on failure
% Pure linear script (no local functions): the MCP wrapper evals script
% text in its own context, which breaks local-function layouts.
% Root detection: the MCP wrapper cd's into the script's own directory
% (tests/) and mfilename('fullFile') may return a bare name, so prefer
% pwd-relative detection.
mf = mfilename('fullFile');
if ~isempty(mf) && any(mf == filesep | mf == '/')
    root = fileparts(fileparts(mf));
else
    % Headless wrapper: cwd is the script's own directory (tests/).
    % GOTCHA (M6): the MCP wrapper's DEFAULT cwd is the scratch sandbox
    % (C:\Users\...\hermes\cache\scratch\matlab-sandbox), which holds a
    % STALE snapshot of tests/ + golden/ - the root detection below then
    % silently globs the snapshot (fewer classes, no error). Always invoke
    % as  cd('D:/Source/hermes-dir/fdedge/tests'); run_tests;  so root
    % lands on the real tree (verified: 67 tests, M6).
    if exist(fullfile(pwd, 'run_tests.m'), 'file')
        root = fileparts(pwd);
    elseif exist(fullfile(pwd, 'tests', 'run_tests.m'), 'file')
        root = pwd;
    elseif exist(fullfile(pwd, 'fdedge', 'tests', 'run_tests.m'), 'file')
        root = fullfile(pwd, 'fdedge');
    else
        error('run_tests:noRoot', 'cannot locate the fdedge project root from cwd %s', pwd);
    end
end
% Explicit addpath (NOT genpath): in this R2024a wrapper genpath returns a
% char array, not a cell array, so addpath(genpath(...)) is a no-op. The
% toolbox is a PLAIN directory (fdedge/), not a +package: headless batch
% mode in this install never resolves package-qualified names (12-line
% minimal repro; see the architecture note in the v1 plan). Add each dir.
addpath(root);
addpath(fullfile(root, 'tests'));
tb = fullfile(root, 'fdedge');
addpath(tb);
subdirs = {'fd', 'est', 'calib', 'tracker', 'report', 'signal'};
for k = 1:numel(subdirs)
    d = fullfile(tb, subdirs{k});
    if exist(d, 'dir')
        addpath(d);
    end
end
if exist(fullfile(root, 'examples'), 'dir')
    addpath(fullfile(root, 'examples'));
end
if exist(fullfile(root, 'golden'), 'dir')
    addpath(fullfile(root, 'golden'));
end
cd(root);
% NOTE: index loop (NOT `for f = dir(...)`) — the headless wrapper
% mis-handles element-wise for-loops over dir() struct arrays (runs the
% body once). Index loops are verified correct.
classes = {};
dd = dir(fullfile(root, 'tests', 'test_*.m'));
for k = 1:numel(dd)
    classes{end+1} = strrep(dd(k).name, '.m', '');
end
nPass = 0; nFail = 0;
for i = 1:numel(classes)
    try
        res = runtests(classes{i});
        nPass = nPass + sum([res.Passed]);
        nFail = nFail + sum([res.Failed]);
        fprintf('%s: %d passed, %d failed\n', classes{i}, sum([res.Passed]), sum([res.Failed]));
    catch ME
        nFail = nFail + 1;
        fprintf('%s: CRASH %s\n', classes{i}, ME.identifier);
        fprintf('   at %s line %d\n', ME.stack(1).name, ME.stack(1).line);
    end
end
fprintf('TOTAL: %d passed, %d failed\n', nPass, nFail);
if nFail > 0
    disp('TESTS_FAILED');
    exit(1);
else
    disp('TESTS_PASSED');
end
