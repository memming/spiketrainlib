% Add all path required to run the library
p = mfilename('fullpath');
[~, idx] = find(p == '/', 1, 'last');
base = p(1:idx);

% add main components
addpath([base 'kernel']);
addpath([base 'generate']);
addpath([base 'util']);
addpath([base 'visual']);
addpath([base 'divergence']);
addpath([base 'metric']);
addpath([base 'regression']);
