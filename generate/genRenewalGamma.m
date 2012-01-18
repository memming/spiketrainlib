function spikeTrains = genRenewalGamma(N, M, param)
% Generate renewal process with interval distribution that follows Gamma dist.
% spikeTrains = genRenewalGamma(N, M, param)
%
% Input
%   N: trials per spikeTrains
%   M: number of sets of trials
%   param.rate: mean rate
%   param.shape: shape parameter
%   param.T: length of simulation
%
% samples from gamma distribution: gamrnd(shape_param,scale_param)
% mean = shape_param * scale_param
% var = shape_param * scale_param^2
% The mean firing rate is 1/mean of interval distribution, rate = 1/mean
%
% $Id: genRenewalGamma.m 32 2011-01-24 20:45:17Z memming $
% Copyright 2010 iocane project. All rights reserved.

% Redistribution and use in source and binary forms, with or without
% modification, are permitted provided that the following conditions are met:
%  - Redistributions of source code must retain the above copyright notice,
%    this list of conditions and the following disclaimer.
%  - Redistributions in binary form must reproduce the above copyright notice,
%    this list of conditions and the following disclaimer in the documentation
%    and/or other materials provided with the distribution.
%  - Neither the name of the iocane project nor the names of its contributors
%    may be used to endorse or promote products derived from this software
%    without specific prior written permission.
%
% THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
% AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
% IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
% ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE
% LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
% CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
% SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
% INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
% CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
% ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
% POSSIBILITY OF SUCH DAMAGE.

spikeTrainsTemplate.duration = param.T;
spikeTrainsTemplate.source = '$Id: genRenewalGamma.m 32 2011-01-24 20:45:17Z memming $';
spikeTrainsTemplate.samplingRate = Inf;
spikeTrainsTemplate.N = N;
spikeTrainsTemplate.data = cell(N, 1);

a = param.shape;
b = 1 / param.rate / a;
meanI = a*b;
stdI = sqrt(a)*b;
nSpikes_mean = (param.T / b);
l95q_interval = gaminv(0.05, a, b); % lower quantile
Nmax = ceil(param.T / l95q_interval);
transientT = 5 * nSpikes_mean; % to make it stationary, simulate a bit more
TT = param.T + transientT;

for kM = 1:M
    spikeTrains(kM) = spikeTrainsTemplate;
    for k = 1:N
	intervals = gamrnd(a, b, Nmax, 1);
	st = cumsum(intervals);
       	tLen = st(end);
	while tLen < TT % If it was not enough, we generate some more
	    Nmax = Nmax * 2;
	    intervals = gamrnd(a, b, Nmax, 1);
	    st = cumsum(intervals);
	    tLen = st(end);
	end
	st = st(st >= transientT & st < TT) - transientT;
	spikeTrains(kM).data{k} = st;
    end
end
