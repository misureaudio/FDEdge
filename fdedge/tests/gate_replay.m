function [pass, detail] = gate_replay(replayFile, outFile, goldenFile, excFile)
% Run a replay driver in a temp dir and gate its output file against the
% golden via gate_files (two-tier gate).
%   replayFile : path to the replay driver script (examples/replay_*.m)
%   outFile    : the results filename the driver writes (cwd-relative)
%   goldenFile : path to the golden snapshot
%   excFile    : exception log for stage-2-only failures (may be '')
if nargin < 4 || isempty(excFile)
    excFile = '';
end
orig = pwd;
tmp = fullfile(tempdir, ['replay_' matlab.lang.makeValidName(outFile)]);
if exist(tmp, 'dir'); rmdir(tmp, 's'); end
mkdir(tmp);
% Copy the driver into tmp and run it from there: this MATLAB wrapper
% executes scripts with cwd = the script's directory, so the driver's
% cwd-relative output file lands in tmp (in stock MATLAB, cd(tmp) below
% already makes cwd = tmp). The driver must therefore be self-contained
% (its own addpath block) - the caller's path is not inherited.
[~, base, ext] = fileparts(replayFile);
base = [base ext];
copyfile(replayFile, fullfile(tmp, base));
cd(tmp);
try
    run(fullfile(tmp, base));
catch ME
    cd(orig);
    pass = false;
    detail = sprintf('replay crashed: %s (line %d in %s)', ME.message, ME.stack(1).line, ME.stack(1).name);
    return
end
cd(orig);
[pass, detail] = gate_files(fullfile(tmp, outFile), goldenFile, excFile);
end
