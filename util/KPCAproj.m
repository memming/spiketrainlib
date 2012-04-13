function [proj2] = KPCAproj(KM, KS, eigVec)
% [proj2] = KPCAproj(KM, KS, eigVec)
% Kernel PCA projection of KS on KM
% Input
%   KM: (NxN) kernel matrix of training set
%   KS: (MxN) gram matrix of test set on training set
%   eigVec: (NxN) result of KPCA.m
%
% Output
%   proj2: (MxN) projection of test set on training set
%
% See also: KPCA

N = size(KM, 1);
M = size(KS, 1);
assert(N == size(KM,2), 'KM should be square');
assert(N == size(KS,2), 'KS should be M x N');

Kc2 = KS - repmat(sum(KM, 1), M, 1) / N - repmat(sum(KS, 2), 1, N) / N ...
    + sum(KM(:)) / N^2;

proj2 = Kc2 * eigVec;
