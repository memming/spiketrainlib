function xidxs = xvalidationLPO(idx1, idx2)
% xidxs = xvalidationLPO(idx1, idx2)
% Leave pair out cross-validation (LPOCV) scheme (Airola et. al. MLSB 2009)
% It is an unbiased way of computing AUC (Cortes et. al. 2007)
%
% Input
%   idx1: (N1x1) vector of indices for class 1
%   idx2: (N1x1) vector of indices for class 2
%
% Output
%   xidxs: {kxValidation x 2} training, test indices
%
% See also: xvalidationIdx
%
% Copyright 2011 Memming. All rights reserved.

nPairs = length(idx1) * length(idx2);
xidxs = cell(nPairs, 2);

k = 0;
for k1 = 1:length(idx1)
    for k2 = 1:length(idx2)
	k = k + 1;
	xidxs{k,1} = [idx1(1:(k1-1)); idx1((k1+1):end); idx2(1:(k2-1)); idx2((k2+1):end)];
	xidxs{k,2} = [idx1(k1); idx2(k2)];
    end
end
