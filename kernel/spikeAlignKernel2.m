function val = spikeAlignKernel2(ks, st1, st2, param)
% Computes the Gram matrix given spike trains st1 and st2, both cells of arrays
% This code is based on interspike interval, the first interval is chosen
%   with respect to firstSpike
% param is the kernel size for the exponential kernel i.e.
%   exp(-|t_1-t_2|*param)
% Adding/Deleting spike is equivalent to setting interval to 0
%   
% Author Sohan Seth (sohan@cnel.ufl.edu)
% Implemented dynamic programming, and adapted to stkernel library by Memming (2011/05/11)
%
% $Id: spikeAlignKernel2.m 78 2011-07-28 21:15:55Z memming $

st1 = st1(:); st2 = st2(:);
% convert to ISI
st1 = st1 - [0; st1(1:end-1)];
st2 = st2 - [0; st2(1:end-1)];
st1 = param(1) * st1;
st2 = param(1) * st2;

N1 = length(st1);
N2 = length(st2);
if N1 == 0 && N2 == 0
    val = 1; return;
elseif N1 == 0 && N2 ~= 0 % Only one is empty ...
    val = exp(-sum(st2)); return;
elseif N2 == 0 && N1 ~= 0 % Only one is empty ...
    val = exp(-sum(st1)); return;
end

% Recursive implementation of alignment kernel 
% Considers mean of all possible paths rather than only the minimum valued path
% Therefore, K(x,x) depends on x, but symmetric K(x,y) = K(y,x) 
%           and *hopefully* strictly positive definite

%% Original implementation by Sohan
% val = mean([spikeAlignKernel2(ks,st1(1:end-1),st2(1:end),param) * kernel(st1(end),[],param), ...
%     spikeAlignKernel2(ks,st1(1:end-1),st2(1:end-1),param) * kernel(st1(end),st2(end),param), ...
%     spikeAlignKernel2(ks,st1(1:end),st2(1:end-1),param) * kernel([],st2(end),param)]);

% Dynamic programming (by Memming)
val = nan(N1+1, N2+1); % val(k1, k2) = spikeAlignKernel(st1(1:k1-1),st2(1:k2-1))
% initialize the boundary values
val(1,1) = 1;
val(1,2:N2+1) = exp(-cumsum(st2));
val(2:N1+1,1) = exp(-cumsum(st1));
% fill in the interior
totalLen = 2;
while true
    for k1 = 1:N1
	if k1 > totalLen + 1
	    break;
	end
	k2 = totalLen - k1;
	if k2 <= 0 || k2 > N2
	    continue;
	end

	v1 = val(k1  ,k2+1) * exp(-st1(k1));
	v2 = val(k1  ,k2  ) * exp(-abs(st1(k1)-st2(k2)));
	v3 = val(k1+1,k2  ) * exp(-st2(k2));
	val(k1+1, k2+1) = (v1 + v2 + v3) / 3;
    end

    if totalLen == N1 + N2 + 1
	break;
    end
    totalLen = totalLen + 1;
end

% assert(~any(isnan(val(:))), 'Values not full!')
val = val(end,end);
