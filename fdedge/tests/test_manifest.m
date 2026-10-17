classdef test_manifest < matlab.unittest.TestCase
    % M11: all four replays reproduce their golden snapshots. For every
    % entry in tests/golden_manifest.txt (file, len, sum of uint8): run
    % the corresponding examples/replay_<profile>.m in a temp dir, check
    % numel/sum against the manifest, and pass the two-tier gate
    % (stage 1 numeric, stage 2 byte-identical). Prints one len/sum line
    % per entry. run_tests cd's to the project root before running the
    % suite.
    methods (Test)
        function allFourReplaysMatchManifest(tc)
            % NOTE: the TestCase handle is 'tc', not 't': run() executes
            % the replay drivers in THIS workspace, and the PoC-derived
            % drivers use 't' as their time-grid variable - naming the
            % handle 't' would let the driver clobber it.
            root = pwd;
            exc = fullfile(root, 'tests', 'replay_exceptions.md');
            % ---- parse the manifest (line by line; skip '#' comments) ----
            txt = fileread(fullfile(root, 'tests', 'golden_manifest.txt'));
            lines = splitlines(txt);
            nman = 0;
            manFile = {}; manLen = zeros(0,1); manSum = zeros(0,1);
            for k = 1:numel(lines)
                L = strtrim(lines{k});
                if isempty(L) || L(1) == '#'
                    continue;
                end
                nman = nman + 1;
                parts = split(strtrim(L));
                manFile{nman} = parts{1};
                manLen(nman) = str2double(parts{2});
                manSum(nman) = str2double(parts{3});
            end
            tc.verifyEqual(nman, 4, 'expected 4 manifest entries');
            % ---- per entry: replay in a temp dir, check, gate ----
            for i = 1:nman
                % golden/results_<profile>.txt -> replay_<profile>.m
                prof = strrep(manFile{i}, 'golden/results_', '');
                prof = strrep(prof, '.txt', '');
                driver = ['replay_' prof '.m'];
                out = ['results_' prof '.txt'];
                golden = fullfile(root, manFile{i});
                tmp = fullfile(tempdir, ['manifest_' driver]);
                if exist(tmp, 'dir')
                    rmdir(tmp, 's');
                end
                mkdir(tmp);
                copyfile(fullfile(root, 'examples', driver), tmp);
                orig = pwd;
                crashed = false;
                errmsg = '';
                try
                    cd(tmp);
                    run(fullfile(tmp, driver));
                catch ME
                    crashed = true;
                    errmsg = ME.message;
                end
                cd(orig);
                if crashed
                    tc.verifyFail(sprintf('%s: replay crashed: %s', manFile{i}, errmsg));
                end
                % ---- len/sum check against the manifest ----
                % The goldens are authoritative snapshots with mixed
                % line-ending conventions (poc: CRLF, the rest: LF); the
                % replay always writes LF. Normalize the replay to the
                % golden's convention before comparing len/sum.
                gfid = fopen(golden, 'r');
                G = fread(gfid, '*uint8');
                fclose(gfid);
                goldenCRLF = any(G == 13);
                fid = fopen(fullfile(tmp, out), 'r');
                b = fread(fid, '*uint8');
                fclose(fid);
                b = b(:)';
                s = char(b);
                crlf = char([13 10]);
                lf = char([10]);
                if goldenCRLF
                    s = strrep(strrep(s, crlf, lf), lf, crlf);
                else
                    s = strrep(s, crlf, lf);
                end
                b = uint8(s);
                gotLen = numel(b);
                gotSum = sum(b, 2, 'double');
                fprintf('manifest: %s  len = %d (want %d)  sum = %d (want %d)\n', ...
                    manFile{i}, gotLen, manLen(i), gotSum, manSum(i));
                tc.verifyEqual(gotLen, manLen(i), ...
                    sprintf('%s: len %d != %d', manFile{i}, gotLen, manLen(i)));
                tc.verifyEqual(gotSum, manSum(i), ...
                    sprintf('%s: sum %d != %d', manFile{i}, gotSum, manSum(i)));
                % ---- two-tier gate ----
                [pass, detail] = gate_files(fullfile(tmp, out), golden, exc);
                tc.verifyTrue(pass, detail);
            end
        end
    end
end
