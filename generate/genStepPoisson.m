function spikeTrains = genStepPoisson(N, M, param)
% Hypothesis test example experiment - Poisson process with chaning rate
% spikeTrains = genStepPoisson(N, M, param)
%
% Input
%   N: trials per spikeTrains
%   M: number of sets of trials
%   param.lambda1, lambda2: mean number of total action potentials
%	(irrespective of duration)
%   param.sectionLength: duration of each homogeneous section
%
%    +-----+                           +-----+
%    |     +-----+         vs.   +-----+     |
%  --+-----+-----+-----> t     --+-----+-----+-----> t
%    lambda1  2                    2      1
%
% $Id: genStepPoisson.m 32 2011-01-24 20:45:17Z memming $
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

sectionLength = param.sectionLength;
T = 2 * sectionLength;
lambda1 = param.lambda1;
lambda2 = param.lambda2;

spikeTrainsTemplate.duration = T;
spikeTrainsTemplate.source = '$Id: genStepPoisson.m 32 2011-01-24 20:45:17Z memming $';
spikeTrainsTemplate.samplingRate = Inf;
spikeTrainsTemplate.N = N;
spikeTrainsTemplate.data = cell(N, 1);

for kM = 1:M
    spikeTrains(kM) = spikeTrainsTemplate;
    for k = 1:N
	n1 = poissrnd(lambda1);
	n2 = poissrnd(lambda2);
	st = sort([rand(n1, 1) * sectionLength; rand(n2, 1) * sectionLength + sectionLength]);
	spikeTrains(kM).data{k} = st;
    end
end

% vim:ts=8:sts=4:sw=4
