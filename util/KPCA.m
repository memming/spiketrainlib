function [proj, eigVal, eigVec] = KPCA(KM, kLargest)
% [proj, eigVal, eigVec] = KPCA(KM)
% Kernel PCA with empirical centering in the feature space
% Input
%   KM: (NxN) kernel matrix of training set
%   kLargest: (1/optional) only obtain the largest eigenvalues
%
% Output
%   proj: (MxN) projection of test set on training set
%   eigVec: (Nx1) eigenvalues
%   eigVec: (NxN) eigenvectors (note that they are not normalized)
%
% See also: KPCAproj

if nargin < 2
    kLargest = min(size(KM,1), 10);
end

N = size(KM, 1);
H = eye(N) - repmat(1/N, N, N);
Kc = H * KM * H; % centered kernel matrix
[eigVec eigVal] = eigs(Kc, kLargest);
eigVal = diag(eigVal);
[eigVal, sidx] = sort(eigVal, 'descend');
eigVec = eigVec(:, sidx);
proj = Kc * eigVec;