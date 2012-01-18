function [candidates] = autoKernelSize2(intensityHandle, aggregationHandle, spikeTrains)
% Suggests kernel sizes for two stage linear scaling kernels.
% Kernel sizes are assumed to be strictly positive.
% These kernels has a cascade of linear and nonlinear operations
%
%  nonlinear3 <- linear2 <- nonlinear2 <- linear1 <- nonlinear1 (st1, st2)
%
% The two stages where the kernel sizes are required are the linear steps.
% The nonlinear1 and nonlinear2 needs to be provided to this function.
%
% kernel = nonLinear3(aggregationHandle(?intensityHandle(st{1}, st{2})/ks1)/ks2)
% TODO - aggregationHandle is separate
%
% Input:
%   intensityHandle: function handle for computing the quantity that needs
%			to be scaled by a linear kernel size
%   aggregationHandle: the result of intensityHandle will be rescaled through
%			
%   spikeTrains: {Nx1} set of spike trains that the kernel will be applied for
% Output:
%   candidates: {9x2} suggested kernel sizes for intensity and aggregation
%
% Copyright 2010-2012 Memming. All rights reserved.

%% Compute the candidate kernel sizes for the INNER function handle
candidates1 = autoKernelSizeScalar(intensityHandle, spikeTrains);
if ~iscell(candidates1)
    candidates = 0; % TODO Change 0 to NaN
    return;
end

%   candidates1: (3x1) suggested kernel sizes for intensity handle
%   candidates2: (3x3) suggested kernel sizes for aggregation handle 
%	    for each intensity handle (each row is for one kernel size 
%	    for intensity handle)
candidates = [];

%% Compute the candidate kernel sizes for the OUTTER function handle
% This is done for each candidate kernel size chosen in the first part
for cIdx = 1:length(candidates1)
    c1 = candidates1{cIdx};
    candidates2 = autoKernelSizeScalar(...
	@(st1, st2)(aggregationHandle(st1, st2, c1)), spikeTrains);
    if ~iscell(candidates2)
	candidates = 0;
	return;
    end
    candidates = ...
	[candidates; [repmat(c1, length(candidates2), 1), cell2mat(candidates2)]];
end

%candidates = sortrows(unique(candidates, 'rows'));
candidates = sortrows(candidates);

candidates = mat2cell(candidates, ones(size(candidates,1),1), size(candidates,2));
