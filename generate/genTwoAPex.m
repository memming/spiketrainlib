function spikeTrains = genTwoAPex(N, M, param)
% Generate two point processes with less than two APs
% spikeTrains = genTwoAPex(N, M, param)
%
% Input
%   N: trials per spikeTrains
%   M: number of sets of trials
%   param.jitter: standard deviation of spike timing jitter
%   param.type: 'correlated' or 'uncorrelated'
%
% PP1. ISI maintained (t1, t1+alpha), t1 ~ unif(0.5,1)
% PP2. (t1, t2), t1 ~ unif(0.5,1), t2 ~ unif(0.5, 1) + alpha
% Both t1 and t2 can be lost with probability p
%
% $Id: genTwoAPex.m 32 2011-01-24 20:45:17Z memming $
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

switch(lower(param.type))
case {'correlated'}
    isCorrelated = 1;
case {'uncorrelated', 'ptst'}
    isCorrelated = 0;
otherwise
    error('Unknown type: use either correlated or renewal');
end

alpha = 0.3;
tWidth = 0.1;
tOffset = 0.2;
jitterSigma = param.jitter; % jitterSigma = 0.01;
duration = 2 * tOffset + tWidth + alpha;
p = 0.1;

lossyAPs = @(st,p)(st(rand(size(st)) >= p));

for kM = 1:M
    spikeTrains(kM).N = N;
    spikeTrains(kM).duration = duration;
    spikeTrains(kM).source = '$Id: genTwoAPex.m 32 2011-01-24 20:45:17Z memming $';
    spikeTrains(kM).data = cell(N, 1);
    spikeTrains(kM).samplingRate = Inf;

    for k = 1:N
	t1a = rand(N, 1) * tWidth + tOffset;
	if ~isCorrelated
	    t1b = rand(N, 1) * tWidth + tOffset;
	    t2b = rand(N, 1) * tWidth + tOffset + alpha + randn(N, 1) * jitterSigma;
	    spikeTrains(kM).data{k} = lossyAPs([t1b(k); t2b(k)], p);
	else
	    t2a = t1a + alpha + randn(N, 1) * jitterSigma;
	    spikeTrains(kM).data{k} = lossyAPs([t1a(k); t2a(k)], p);
	end
    end
end

% vim:ts=8:sts=4:sw=4
