function [proj, eigVal, eigVec] = KPCA(KM)
% [proj, eigVal, eigVec] = KPCA(KM)
% Kernel PCA with empirical centering in the feature space
% Input
%   KM: (NxN) kernel matrix of training set
%
% Output
%   proj: (MxN) projection of test set on training set
%   eigVec: (Nx1) eigenvalues
%   eigVec: (NxN) eigenvectors (note that they are not normalized)
%
% See also: KPCAproj

N = size(KM, 1);
H = eye(N) - repmat(1/N, N, N);
Kc = H * KM * H; % centered kernel matrix
[eigVec eigVal] = eig(Kc);
proj = Kc * eigVec;
eigVal = diag(eigVal);
