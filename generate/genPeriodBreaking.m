function [sts phaseIdx periodIdx] = genPeriodBreaking(param, N)
% Generate spike trains with a periodic firing, and additional spike inside
% with a constant breaking ratio....
%
%        period
% |<---------------->|
% |            |     |
% |<---------->|<--->|
%    period*r

k = 1;
sts = cell(length(param.periodList) * length(param.phaseList), 1);
for kPeriod = 1:length(param.periodList)
    for kPhase = 1:length(param.phaseList)
	period = param.periodList(kPeriod);
	n = ceil(param.T / period);
	st = period * param.phaseList(kPhase) + (0:n-1) * period;
	sts{k} = st(st < param.T);
	phaseIdx(k) = kPhase;
	periodIdx(k) = kPeriod;
	k = k + 1;
	if k > N; break; end
    end
    if k > N; break; end
end
sts = sts(1:N);
