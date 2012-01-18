function spikeTrains = genHomogeneousPoisson(N, M, param)
% Generate homogeneous Poisson process for testing
% spikeTrains = genHomogeneousPoisson(N, M, param)
%
% Input
%   N: trials per spikeTrains
%   M: number of sets of trials
%   param.lambda: mean number of total action potentials
%	(irrespective of duration)
%   param.tOffset: offset at the beginning
%   param.duration: duration of spiking
%
% $Id: genHomogeneousPoisson.m 32 2011-01-24 20:45:17Z memming $
% Copyright 2009 Memming. All rights reserved.

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

lambda = param.lambda;
tOffset = param.tOffset;
duration = param.duration;

spikeTrainsTemplate.N = N;
spikeTrainsTemplate.duration = 2 * tOffset + duration;
spikeTrainsTemplate.source = '$Id: genHomogeneousPoisson.m 32 2011-01-24 20:45:17Z memming $';
spikeTrainsTemplate.data = cell(N, 1);
spikeTrainsTemplate.samplingRate = Inf;

for kM = 1:M
    spikeTrains(kM) = spikeTrainsTemplate;

    for k = 1:N
	spikeTrains(kM).data{k} = ...
	    tOffset + sort(rand(poissrnd(lambda), 1)) * duration;
    end
end

% vim:ts=8:sts=4:sw=4
