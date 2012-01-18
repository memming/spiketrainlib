function [subsamples, targets, sts, temporalWindow] = timeSeriesRandomSample(x, st, samplingFreq, nSubsample, subsampleInterval, isRandom, temporalWindow)
% Randomly subsample a time series and associated spike train for regression.
% By random, we mean uniformly random in time.
% A causal window of spike trains are obtained.
%
% Input
%   x: time series sampled at samplingFreq
%   st: spike train in continuous time
%   samplingFreq: sampling frequency
%   nSubsample: number of samples to obtain
%   subsampleInterval: [start_time end_time]
%   temporalWindow: length of the temporal window to grab the spike trains
%		    By default, 3 times the 95% quantile of ISI is chosen.
%
% Output
%   subsamples: times of random samples
%   targets: samples of x obtained at rounded subsample times
%   sts: samples of spike trains obtained in [-temporalWindow 0] + subsamples
%   temporalWindow: length of the temporal window that is used
%
% $Id: timeSeriesRandomSample.m 94 2011-11-03 14:45:16Z memming $

if nargin < 6; isRandom = true; end

if isRandom
    subsamples = sort(rand(nSubsample, 1) * diff(subsampleInterval) ...
        + subsampleInterval(1));
else
    subsamples = (1:nSubsample)'/nSubsample * diff(subsampleInterval) ...
        + subsampleInterval(1);
end

% estimate required temporal window size
if nargin < 7
    temporalWindow = quantile(diff(st), 0.95) * 3; % Also try 5 times :P
end
% chop up spike trains and make a regression data set
sts = cell(nSubsample, 1);
for sidx = 1:nSubsample
    t = subsamples(sidx);
    sts{sidx} = st(st >= t - temporalWindow & st < t) - t + temporalWindow;
    % -t: spike timing starts at zero
    % -temporalWindow: shift spike trains back for temporal window 
end

if ~isempty(x)
    targets = x(round(subsamples * samplingFreq));
else
    targets = [];
end
