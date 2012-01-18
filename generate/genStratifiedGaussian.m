function spikeTrains = genStratifiedGaussian(N, M, param)
% Generate point process for a given stratified specification
% spikeTrains = genStratifiedGaussian(N, M, param)
% The joint distributions are restricted to be multi-variate Gaussian
%
% Input
%   N: trials per spikeTrains
%   M: number of sets of trials
%   param.D: maximum dimension
%   param.T: length of the spike train
%   param.mean: {Dx1} cell each containing d dimensional vector of means
%   param.sigma: {Dx1} cell each containing dxd covaraince matirx
%
% $Id: genStratifiedGaussian.m 32 2011-01-24 20:45:17Z memming $
% Copyright 2010 Memming. All rights reserved.

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

D = param.D;
cumP = cumsum(param.p);
T = param.T;

spikeTrainsTemplate.N = N;
spikeTrainsTemplate.duration = T;
spikeTrainsTemplate.source = '$Id: genStratifiedGaussian.m 32 2011-01-24 20:45:17Z memming $';
spikeTrainsTemplate.data = cell(N, 1);
spikeTrainsTemplate.samplingRate = Inf;

for kM = 1:M
    spikeTrains(kM) = spikeTrainsTemplate;
    for n = 1:N
	d = find(rand < cumP, 1, 'first');
	st = mvnrnd(param.mu{d}, param.sigma{d});
	spikeTrains(kM).data{n} = sort(st);
    end
end

% vim:ts=8:sts=4:sw=4
