function v = stcellfun(fhandle, stc)
% v = stcellfun(fhandle, stc)
%
% Apply the function handle for each cell up to 2 dimensions
% Empty cell content will result in NaN
%
% Copyright 2010-2012 Memming.

v = nan(size(stc));

for k1 = 1:size(stc,1)
    for k2 = 1:size(stc,2)
	if ~isempty(stc{k1,k2})
	    v(k1,k2) = fhandle(stc{k1,k2});
	end
    end
end
