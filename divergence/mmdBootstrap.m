function [rejected, mmd2, threshold, d2] = mmdBootstrap(KM, idx1, idx2, alpha)
% [rejected, mmd2, threshold, d2] = mmdBootstrap(KM, idx1, idx2, alpha)
% Two sample hypothesis test using MMD for a given kernel matrix.
% Note that this function is not specifically designed for spike trains.
% Null hypothesis is that they came from the same distribution.
% Compute MMD and it's associated statistics using bootstrapping.
%
% We bootstrap MMD^2, and test if 0 is in the lower tail.
% (Recommended by Hall and Wilson 1991)
% Resampling is done with replacement (maintain sample size using repeats).
%
% Note that this is uses the biased estimator version of MMD.
%
% Input
%   KM: (NxN) symmetric pd kernel matrix
%   idx1: indices for samples from distribution 1
%   idx2: indices for samples from distribution 2
%   alpha: (1) significance level (default: 0.05)
%
% Output
%   rejected: null hyphothesis rejected? true or false
%   mmd2: statistic MMD^2
%   threshold: the bootstrapped threshold for mmd2
%
% Ref: 
% [1] A. Gretton, K. Borgwardt, M. Rasch, B. Schoelkopf, and A. Smola.
%     A kernel method for the two-sample-problem. NIPS 19, 2007
% [2] A. Gretton, K. Fukumizu, Z. Harchaoui, and B. Sriperumbudur.
%     A fast, consistent kernel two-sample test. NIPS 22, 2009
%
% Copyright 2011 Memming. All rights reserved.

N = size(KM, 1); assert(size(KM,2) == N);

% convert binary index to numerical index
if any(idx1 == 0); idx1 = find(idx1); idx2 = find(idx2); end

N1 = length(idx1); N2 = length(idx2);
assert(N == N1 + N2);
o1 = ones(N1, 1); o2 = ones(N2, 1); o = ones(N, 1);

% MMD^2 = 1/N1^2 \sum_i,j K(x_i, x_j) 
%	+ 1/N1^2 \sum_i,j K(x_i, x_j) 
%	- 2/N1/N2 \sum_i,j K(x_i, y_j) 
mmd2 = o1' * KM(idx1, idx1) * o1 / (N1^2) ...
     + o2' * KM(idx2, idx2) * o2 / (N2^2) ...
     - o1' * KM(idx1, idx2) * o2 / (N1*N2) ...
     - o2' * KM(idx1, idx2)' * o1 / (N1*N2);

%% number of bootstrap samples to generate
% @TechReport{RePEc:qed:wpaper:1127,
%   author={James G. MacKinnon},
%   title={Bootstrap Hypothesis Testing},
%   year=2007,
%   month=Jun,
%   institution={Queen's University, Department of Economics},
%   type={Working Papers},
%   url={http://ideas.repec.org/p/qed/wpaper/1127.html},
%   number={1127},
%   abstract={},
%   keywords={}
% }
B = 9999;
assert(rem(alpha * (B+1),1) == 0, 'This must be an integer');

bmmd2 = zeros(B, 1);
for kBootstrap = 1:B
    bidx1 = idx1(ceil(rand(N1, 1) * N1));
    bidx2 = idx2(ceil(rand(N2, 1) * N2));
    bidx = [bidx1; bidx2];
    pidx = randperm(N);
    bidx = bidx(pidx);
    bidx1 = bidx(1:N1);
    bidx2 = bidx((N1+1):end);

    bmmd2(kBootstrap) = ...
	  o1' * KM(bidx1, bidx1) * o1 / (N1^2) ...
	+ o2' * KM(bidx2, bidx2) * o2 / (N2^2) ...
	- o1' * KM(bidx1, bidx2) * o2 / (N1*N2) * 2;
end

% Compute one of the following, and take the alpha quantile
% (a) (\hat{\theta}^\ast - \hat{\theta})
% (b) (\hat{\theta}^\ast - \hat{\theta}) / \hat{\sigma}^\ast    (Pivoting)
% d2 = d2 / std(d2); % pivoting is wrong, because I don't know the std of \hat{theta} itself, I cannot normalize it.

d2 = bmmd2;
threshold = quantile(d2, 1 - alpha);

if mmd2 > threshold
    rejected = true;
else
    rejected = false;
end
