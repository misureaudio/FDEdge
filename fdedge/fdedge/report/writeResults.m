function writeResults(fid, lines)
%WRITERESULTS Write preformatted result lines to an open file.
%   writeResults(fid, lines) writes each cell element of LINES to the
%   file behind file ID FID, verbatim, followed by one newline.
%
% The replay drivers keep their verbatim PoC fprintf blocks (the byte
% gate requires it); this helper serves new consumers such as
% examples/demo.m, which format their blocks with sprintf and hand them
% here in one call.
for k = 1:numel(lines)
    fprintf(fid, '%s\n', lines{k});
end
end
