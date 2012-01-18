function [pvalues, results, resultsTC, mmd2s] = divMMDsupBootstrap(sts1, sts2, kernelList, T, alpha, verbose)
% [pvalues, results, resultsTC, mmd2s] = divMMDsupBootstrap(sts1, sts2, kernelList, T, alpha, verbose)
% sup-Maximum Mean Discrepancy test (Sriperumbudur et al., 2009).
% Runs MMD for different kernel sizes (using autoParam), then takes the
% maximum as the statistic for classical hypothesis testing.
% When multiple kernels are provided, each kernel is treated separately.
%
% Input
%   sts1, sts2: cell array of spike trains
%   kernelList: list of kernels (kernelStruct from kernelFactory or string)
%   T: length of spike trains
%   alpha: size of test (significance level of p-value)
%
% Output
%   pvalues: estimated p-value via bootstrap
%   results: rejection of null or not
%   resultsTC: computation time
%   mmd2s: the test statistic
%
% Ref: Bharath Sriperumbudur, Kenji Fukumizu, Arthur Gretton, Gert Lanckriet,
% Bernhard Schoelkopf. Kernel choice and classifiability for RKHS embeddings
% of probability distributions. In Advances in Neural Information Processing
% Systems 22, pages 1750–1758. 2009.
%
% Copyright 2011 Memming. All rights reserved.

if nargin < 5 || isempty(alpha); alpha = 0.05; end
if nargin < 6 || isempty(verbose); verbose = false; end

assert(length(sts1) == length(sts2));
idx1 = 1:length(sts1);
idx2 = length(sts1)+(1:length(sts2));
sts = cell(length(sts1) + length(sts2), 1);
sts(idx1) = sts1;
sts(idx2) = sts2;

%% Setup kernel parameters
rs = rand('seed'); % <-- !!!

for kKernel = 1:length(kernelList)
    if isstr(kernelList{kKernel}) && strcmp(kernelList{kKernel}, 'wilcoxon')
	tic;
	totalCount = cellfun('length', sts);
	[pvalues(kKernel) rejected stats] = ...
		ranksum(totalCount(idx1), totalCount(idx2), 'alpha', alpha);
	results(kKernel) = rejected;
	resultsTC(kKernel) = toc;
	mmd2s(kKernel) = stats.ranksum; % NOT really MMD2
	if verbose; 
	    fprintf('%20s [%d]     p-value [%f]\n', 'Wilcoxon ranksum', 0, pvalues(kKernel));
       	end
	continue;
    end

    if isstr(kernelList{kKernel})
	ks = kernelFactory(kernelList{kKernel}, T);
    else
	ks = kernelList{kKernel};
    end
    ksizeList = ks.autoParam(ks, sts);
    if ~iscell(ksizeList)
	error('kernel size list is not a cell!');
    end
    nks = length(ksizeList);
    if verbose; fprintf('%20s [%d] ', ks.name, nks); end
    
    %% Compute kernel matrix
    mmd2All = zeros(nks, 1);
    d2All = zeros(nks, 9999); % TODO magic number in mmdBootstrap
    tic;
    for ksidx = 1:size(ksizeList,1)
	KM = computeKernelMatrix(ks, sts, ksizeList{ksidx});

	%% Perform hypothesis test
	rand('seed', rs); % hack to get the same bootstrap samples (TODO)
	[rejected, mmd2_instance, threshold, d2_instance] = ...
	    mmdBootstrap(KM, idx1, idx2, alpha);
	mmd2All(ksidx) = mmd2_instance;
	d2All(ksidx, :) = d2_instance;
	if verbose && rejected
	    fprintf('.');
	elseif verbose
	    fprintf('o');
	end
    end % ksize
    tCost = toc;

    %% Compute the sup MMD value
    mmd2 = max(mmd2All); d2 = max(d2All, [], 1);

    % I want the p-value at which the test is rejected
    p1 = sum(d2 > mmd2) / length(d2);
    p2 = 1 - sum(d2 < mmd2) / length(d2);
    p = (p1 + p2) / 2;

    if verbose; fprintf('   p-value [%f]\n', p); end
    
    pvalues(kKernel) = p;
    results(kKernel) = (p <= alpha);
    resultsTC(kKernel) = tCost;
    mmd2s(kKernel) = mmd2;
end % kKernelList
