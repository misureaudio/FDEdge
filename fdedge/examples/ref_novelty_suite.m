%% REF_NOVELTY_SUITE — reference cases for the spectral-novelty co-gate
%
% Synthetic reference suite for est/spectralNovelty.m (the co-gate
% primitive, no position scan). M = 1024, Fs = 48 kHz, edge at iE = 512
% (window center). Expected: stationary cases ~ 0, edge cases >> 0.1.

root = 'D:/Source/hermes-dir/fdedge';
addpath(fullfile(root, 'fdedge'));
addpath(fullfile(root, 'fdedge', 'fd'));
addpath(fullfile(root, 'fdedge', 'est'));
addpath(fullfile(root, 'fdedge', 'calib'));
addpath(fullfile(root, 'fdedge', 'tracker'));
addpath(fullfile(root, 'fdedge', 'signal'));

Fs = 48000;
M  = 1024;
n  = M + 1;
xs = (0:M)';           % sample indices (1-based content at index k+1)
rng(0);                % deterministic noise for the noisy cases

w = 0.5*(xs >= 512) + 0.02*randn(n,1);
w = w - median(w(1:round(0.25*n)));
[nov,p] = spectralNovelty(w, 512, Fs);
fprintf('real step + noise (sig 0.5, n 0.02):  nov = %8.4f  (p=%d)\n', nov, p);

w = 0.5*(xs >= 512);
w = w - median(w(1:round(0.25*n)));
[nov,p] = spectralNovelty(w, 512, Fs);
fprintf('clean step (0 -> 0.5):                nov = %8.4f  (p=%d)\n', nov, p);

w = 0.2*cos(2*pi*60*xs/Fs - 0.7);
w = w - median(w(1:round(0.25*n)));
[nov,p] = spectralNovelty(w, 512, Fs);
fprintf('60Hz tone (stationary):               nov = %8.4f  (p=%d)\n', nov, p);

w = 0.1*cos(2*pi*1000*xs/Fs + 0.3);
w = w - median(w(1:round(0.25*n)));
[nov,p] = spectralNovelty(w, 512, Fs);
fprintf('1kHz tone (stationary):               nov = %8.4f  (p=%d)\n', nov, p);

w = 0.2*cos(2*pi*60*xs/Fs - 0.7) + 0.5*(xs >= 512);
w = w - median(w(1:round(0.25*n)));
[nov,p] = spectralNovelty(w, 512, Fs);
fprintf('60Hz tone + step at edge:             nov = %8.4f  (p=%d)\n', nov, p);

w = 0.2*cos(2*pi*60*xs/Fs - 0.7) + 0.15*cos(2*pi*400*xs/Fs).*(xs >= 512);
w = w - median(w(1:round(0.25*n)));
[nov,p] = spectralNovelty(w, 512, Fs);
fprintf('60Hz tone + new 400Hz at edge:        nov = %8.4f  (p=%d)  [DEGENERATE: pure single-tone background, see docs]\n', nov, p);

w = 0.2*cos(2*pi*60*xs/Fs - 0.7) + 0.01*randn(n,1) + 0.15*cos(2*pi*400*xs/Fs).*(xs >= 512);
w = w - median(w(1:round(0.25*n)));
[nov,p] = spectralNovelty(w, 512, Fs);
fprintf('60Hz tone + noise(0.01) + new 400Hz:  nov = %8.4f  (p=%d)  [noise breaks the degeneracy]\n', nov, p);

w = 0.2*cos(2*pi*60*xs/Fs - 0.7);
w(513) = w(513) + 0.3;              % single-sample click at the first post sample
w = w - median(w(1:round(0.25*n)));
[nov,p] = spectralNovelty(w, 512, Fs);
fprintf('60Hz tone + 1-sample click at edge:   nov = %8.4f  (p=%d)\n', nov, p);

w = [zeros(512,1); 0.5*ones(513,1)];
[nov,p] = spectralNovelty(w, 512, Fs);
fprintf('silence -> step:                      nov = %8.4f  (p=%d)\n', nov, p);

fprintf('REF_SUITE_DONE\n');
