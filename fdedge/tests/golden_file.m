function g = golden_file(name)
% golden_file(name) - full path to a file in the golden/ directory.
% NOTE: mfilename('fullFile') is unreliable in this headless environment
% (may return a bare name), so anchor on the current working directory
% first: the test runner cd's to the project root before runtests.
cands = { fullfile(pwd, 'golden', name), ...
          fullfile(pwd, 'fdedge', 'golden', name) };
mf = mfilename('fullFile');
if ~isempty(mf) && any(mf == filesep | mf == '/')
    cands{end+1} = fullfile(fileparts(mf), '..', 'golden', name);
end
g = '';
for i = 1:numel(cands)
    if exist(cands{i}, 'file')
        g = cands{i};
        return
    end
end
error('golden_file:notFound', 'golden file %s not found (tried: %s)', name, strjoin(cands, ' ; '));
end
