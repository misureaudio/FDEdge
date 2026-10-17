function [pass, detail] = gate_files(replayFile, goldenFile, excFile)
% TWO-TIER REPLAY GATE on two existing files (v1, per reviews r0/r1).
%   Stage 1 (PRIMARY, portable): every printed number in the replay
%     output agrees with the golden within 2x the place value of its
%     last printed digit (printed integers: exact).
%   Stage 2 (SECONDARY, this environment): byte-identical files.
%     If stage 1 passes but stage 2 fails, triage: strip CR from both
%     and re-compare (pure line-ending difference -> environmental);
%     otherwise log the first difference in excFile for review at M12.
% pass = stage1 && triage-ok.  A genuine numeric difference always fails
% stage 1 and is red regardless of triage.
if nargin < 3 || isempty(excFile)
    excFile = '';
end
rb = readBytes(replayFile);
gb = readBytes(goldenFile);
byteOK = isequal(rb, gb);
% gb(:) first: readBytes orientation varies by environment (fread returns a
% column in this MATLAB wrapper, a row in stock MATLAB); regexp needs a row.
gN = gate_parse(char(gb(:))');
rN = gate_parse(char(rb(:))');
if numel(gN) ~= numel(rN)
    pass = false;
    detail = sprintf('numeric count mismatch: golden %d vs replay %d', numel(gN), numel(rN));
    return
end
for i = 1:numel(gN)
    if abs(gN(i).v - rN(i).v) > gN(i).tol
        pass = false;
        detail = sprintf('stage1 FAIL at number %d: golden %.12g vs replay %.12g (tol %g)', ...
            i, gN(i).v, rN(i).v, gN(i).tol);
        return
    end
end
if byteOK
    pass = true;
    detail = 'stage1+stage2 PASS (byte-identical)';
    return
end
% Triage: byte difference despite numeric agreement.
gLf = normalizeLf(gb);
rLf = normalizeLf(rb);
if isequal(gLf, rLf)
    pass = true;
    detail = 'stage1 PASS, stage2 FAIL -> triaged: line-endings only (CRLF)';
else
    off = firstDiffOffset(gb, rb);
    pass = true;                                    % stage 1 is the gate; record for review
    if ~isempty(excFile)
        fid = fopen(excFile, 'a');
        if fid > 0
            fprintf(fid, '%s: stage2 byte diff at offset %d (stage1 passed)\n', ...
                fileparts(replayFile), off);
            fclose(fid);
        end
    end
    detail = sprintf('stage1 PASS, stage2 FAIL -> logged at offset %d', off);
end
end

function b = readBytes(p)
fid = fopen(p, 'r');
b = fread(fid, '*uint8')';
fclose(fid);
end

function s = normalizeLf(b)
s = char(b)';
s(s == char(13)) = [];
end

function off = firstDiffOffset(a, b)
if numel(a) ~= numel(b)
    off = 0;
    return
end
off = numel(a) + 1;
for i = 1:numel(a)
    if a(i) ~= b(i)
        off = i;
        break
    end
end
end
