function [sts phaseIdx countIdx] = genPeriodic(param, N)
% Generate periodic spike trains with different phases
%
% <---------- param.T ---------->
% |   |         ...          |
%  |   |        ...           |
%   |   |       ...            |
%    |   |      ...             |
%
% $Id$

k = 1;
assert(rem(param.nPhase, 1) == 0);
assert(param.nPhase * length(param.countList) == N);

sts = cell(length(param.countList) * param.nPhase, 1);
for kCount = 1:length(param.countList)
    count = param.countList(kCount);
    period = param.T / count;
    phaseMax = param.T / count;
    for kPhase = 1:param.nPhase
	phase = phaseMax * (kPhase-1) / param.nPhase;
	sts{k} = (0:count-1) * period + phase;
	phaseIdx(k) = kPhase;
	countIdx(k) = kCount;
	k = k + 1;
	if k > N; break; end
    end
    if k > N; break; end
end
sts = sts(1:N);
