% add various paths needed

p = mfilename('fullpath');
[~, idx] = find(p == '/', 1, 'last');
base = p(1:idx);

% add libraries
addpath([base 'lib/PYMentropy/src/']);

% add CDM
addpath([base 'src/']);
