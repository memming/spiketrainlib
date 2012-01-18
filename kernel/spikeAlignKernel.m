function val = spikeAlignKernel(ks, st1, st2, tau)
% val = spikeAlignKernel(ks, st1, st2, tau)
% similar to VP, e.g. exp(-|t_1-t_2|/tau) is the *gain* of moving a
% spike from t_1 to t_2, whereas the *gain* of adding/deleting is exp(-1)
% We use gain since we are interested in similarity and not dissimilarity.
% tau controls the relative contribution of moving and deleting/adding.
%
% NOTE: This kernel is SUSPECTED to be SPD.
%
% Input
%   ks: kernel structure (not used)
%   st1, st2: spike trains (sorted real valued vector)
%   tau: inverse cost for shifting (time constant)
%
% $Id: spikeAlignKernel.m 54 2011-03-25 13:21:21Z memming $
% Author Sohan Seth (sohan@cnel.ufl.edu)
% Implemented dynamic programming, and adapted to stkernel library by Memming (2011/03/04)

N1 = length(st1);
N2 = length(st2);
if N1 == 0 || N2 == 0
    val = exp(-(N1+N2));
    return
end

% Recursive implementation of alignment kernel 
% Considers mean of all possible paths rather than only the minimum valued path
% Therefore, K(x,x) depends on x, but symmetric K(x,y) = K(y,x) 
%           and *hopefully* strictly positive definite
% val = mean([spikeAlignKernel(ks, st1(1:end-1),st2(1:end),tau) * kernel(st1(end),[],tau), ...
%     spikeAlignKernel(ks, st1(1:end-1),st2(1:end-1),tau) * kernel(st1(end),st2(end),tau), ...
%     spikeAlignKernel(ks, st1(1:end),st2(1:end-1),tau) * kernel([],st2(end),tau)]);

% Dynamic programming (by Memming)
val = nan(N1+1, N2+1); % val(k1, k2) = spikeAlignKernel(st1(1:k1-1),st2(1:k2-1))
% initialize the boundary values
val(1,1) = 1;
val(1,2:N2+1) = exp(-(1:N2));
val(2:N1+1,1) = exp(-(1:N1));
% fill in the interior
totalLen = 2;
while true
    for k1 = 1:(N1+1)
	if k1 - 1 > totalLen
	    break;
	end
	k2 = totalLen - k1;
	if k2 <= 0 || k1 <= 0 || k2 > N2 || k1 > N1
	    continue;
	end

	v1 = val(k1  ,k2+1) * kernel(st1(k1), [], tau);
	v2 = val(k1  ,k2  ) * kernel(st1(k1), st2(k2), tau);
	v3 = val(k1+1,k2  ) * kernel([], st2(k2), tau);
	val(k1+1, k2+1) = (v1 + v2 + v3) / 3;
    end

    if totalLen == N1 + N2 + 1
	break;
    end
    totalLen = totalLen + 1;
end

% assert(~any(isnan(val(:))), 'Values not full!')
val = val(end,end);

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function val = kernel(s,t,tau)
% Strictly positive definite kernel over SxS where S can have either zero
% or one element, e.g. s = [], t = [1] or s = [0], t = [2] or s = [3], t = []

if ~isempty(s) && ~isempty(t)
    val = exp(-abs(s-t) / tau);
else
    % assert(~(isempty(s) && isempty(t)), 'Error in implementation');
    val = exp(-1);
end

end
