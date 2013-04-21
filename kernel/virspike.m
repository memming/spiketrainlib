function val = virspike(ks, s1, s2, param)

% Normalize the spike trains
s1 = s1 / ks.T; s2 = s2 / ks.T;

if any(s1 >= ks.T | s1 <= 0)  || any(s2 >= ks.T | s2 <= 0)
    error('spike train is not within the specified time span');
end

% insert virtual spikes at 0 and 1
if ~isempty(s1)
    t1 = zeros(2*length(s1)+1,1);
    t1(end) = 1; t1(2:2:end-1) = s1;
    t1(3:2:end-1) = 0.5*(s1(1:end-1) + s1(2:end));
    a1 = zeros(length(t1),1);
    a1(3:2:end-1) = 0.5*(s1(2:end) - s1(1:end-1));
    a1(1) = s1(1); a1(end) = 1 - s1(end);
else
    t1(1,1) = 0; t1(2,1) = 1;
    a1(1,1) = 1; a1(2,1) = 1;
end

if ~isempty(s2)
    t2 = zeros(2*length(s2)+1,1);
    t2(end) = 1; t2(2:2:end-1) = s2;
    t2(3:2:end-1) = 0.5*(s2(1:end-1) + s2(2:end));
    a2 = zeros(length(t2),1);
    a2(3:2:end-1) = 0.5*(s2(2:end) - s2(1:end-1));
    a2(1) = s2(1); a2(end) = 1 - s2(end);
else
    t2(1,1) = 0; t2(2,1) = 1;
    a2(1,1) = 1; a2(2,1) = 1;
end

t = sort(unique([t1;t2]));

%% The following is equivalent to
% a = interp1(t1,a1,t);
% b = interp1(t2,a2,t);
a = nan(size(t)); b = nan(size(t));
before1 = 1; before2 = 1;
%beforeArr1 = zeros(size(t));
%for t1idx = 1:numel(t1)
%    idx = t1(t1idx) <= t;
%    beforeArr1(idx) = beforeArr1(idx) + 1;
%end
for tidx = 1:numel(t)
    while t1(before1) <= t(tidx)
	before1 = before1 + 1;
	if before1 == numel(t1) + 1
	    break;
	end
    end
    before1 = before1 - 1;
    %disp([before1 beforeArr1(tidx)])
    %before1 = beforeArr1(tidx);
    %before1 = find(t1 <= t(tidx), 1, 'last');
    if t1(before1) == t(tidx)
	a(tidx) = a1(before1);
    else
	if a1(before1) == 0 % increasing
	    a(tidx) = t(tidx) - t1(before1);
	else % decreasing
	    a(tidx) = a1(before1) - (t(tidx) - t1(before1));
	end
    end

    while t2(before2) <= t(tidx)
	before2 = before2 + 1;
	if before2 == numel(t2) + 1
	    break;
	end
    end
    before2 = before2 - 1;
    %before2 = find(t2 <= t(tidx), 1, 'last');
    if t2(before2) == t(tidx)
	b(tidx) = a2(before2);
    else
	if a2(before2) == 0 % increasing
	    b(tidx) = t(tidx) - t2(before2);
	else % decreasing
	    b(tidx) = a2(before2) - (t(tidx) - t2(before2));
	end
    end
end

%{
figure(16); clf;
ax1 = subplot(3,1,1); line(t1, a1); axis equal;
title('spike train 1');
ax2 = subplot(3,1,2); line(t2, a2); axis equal;
title('spike train 2');
ax3 = subplot(3,1,3); hold all; plot(t, a, 'o-'); plot(t, b, 'x-'); axis equal;
title('intermediate computation');
linkaxes([ax1 ax2 ax3]);
xlim([0 1]);
%}

%%
m1 = (a(2:end) - a(1:end-1)) ./ (t(2:end) - t(1:end-1));
m2 = (b(2:end) - b(1:end-1)) ./ (t(2:end) - t(1:end-1));
b1 = a(1:end-1); b2 = b(1:end-1);

val = sum((m1.*m2).*(t(2:end) - t(1:end-1)).^3/3 ...
        + (m1.*b2 + m2.*b1).*(t(2:end) - t(1:end-1)).^2/2 ...
        + (b1.*b2).*(t(2:end) - t(1:end-1)));

assert(~isnan(val), 'NaN is impossible!')
