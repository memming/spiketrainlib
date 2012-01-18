function [sts, feature] = singleFeaturePointProcessFactory(choice, param)
% Returns two point processes a single interpretable statistic difference
% while keeping all other statistic constant. Using this we want to prove 
% that certain divergence are sensitive to the feature under consideration.
% Note that these features are NOT independent, and there are many ways to
% generate datasets with such features. So the word "single" used in the
% function is not really appropriate, but it represents the spirit of the work.
%
% Input (if no input is given, the total number of choices are returned)
%   choice: an integer that indicates which point process to generate
%   param: parameters for the generation, depends on each choice
%	param.T: duration of the process
%	param.N: how many trials to generate
%	for other parameters see the code
%
% Features
%   'count distribution'
%
%
% $Id: singleFeaturePointProcessFactory.m 102 2011-11-06 21:00:25Z memming $

if nargin < 1
    sts = 3;
end

switch choice
    case 1
	feature = 'variance of count (negative binomial)';
	% Count distribution follows negative binomial with same mean
	% but different variance.
	sts = cell(param.N, 2);
	sts(:, 1) = genNegBinHomogeneous(param.meanCount, param.varCount(1), param.T, param.N);
	sts(:, 2) = genNegBinHomogeneous(param.meanCount, param.varCount(2), param.T, param.N);
    case 2
	feature = 'symmetric step PSTV';
	sts = cell(param.N, 2);
	stsParts = cell(param.N, 2);
	stsParts(:, 1) = genNegBinHomogeneous(param.meanCount, param.varCount(1), param.T/2, param.N);
	stsParts(:, 2) = genNegBinHomogeneous(param.meanCount, param.varCount(2), param.T/2, param.N);
	pidx1 = randperm(param.N); pidx2 = randperm(param.N);
	for k = 1:param.N
	    sts{k, 1} = [stsParts{k, 1}; param.T/2 + stsParts{pidx1(k), 2}];
	    sts{k, 2} = [stsParts{k, 2}; param.T/2 + stsParts{pidx2(k), 1}];
	end
    case 3
	feature = 'mean count with Poisson';
	for k = 1:param.N
	    sts{k, 1} = sort(rand(poissrnd(param.lambda(1)), 1)) * param.T;
	    sts{k, 2} = sort(rand(poissrnd(param.lambda(2)), 1)) * param.T;
	end
    otherwise
	error('Unknown option');
end
