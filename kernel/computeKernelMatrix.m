function KM = computeKernelMatrix(kernelStruct, sts, ksize, sts2)
% KM = computeKernelMatrix(kernelStruct, sts, ksize, sts2)
% Computes the kernel matrix.
%
% Input:
%   kernelStruct: kernel structure (obtained via kernelFactory)
%   sts: cell array of spike trains
%   ksize: kernel size vector
%   sts2: (optional) if given computes KM(k1, k2) = ks.kernel(sts{k1}, sts{k2});
% Output:
%   KM: kernel matrix
%
% For certain kernels, the kernelStructure has a method 'KM' which could be
% faster (automatically called). Otherwise it is O(n^2).
%
% Copyright 2010-2012 Memming

N = length(sts);
if nargin < 4
    sts2 = sts;
    % if special kernel matrix operation is already built in or not
    if isfield(kernelStruct, 'KM') 
	KM = kernelStruct.KM(kernelStruct, sts, ksize);
	return;
    end
    if isfield(kernelStruct, 'isHermitian') && kernelStruct.isHermitian
	KM = zeros(N);
	for k1 = 1:N
	    for k2 = k1:N
		KM(k1, k2) = kernelStruct.kernel(kernelStruct, ...
						sts{k1}, sts{k2}, ksize);
		KM(k2, k1) = conj(KM(k1, k2));
	    end
	end
	return;
    end
end

N2 = length(sts2);
KM = zeros(N, N2);
for k1 = 1:N
    for k2 = 1:N2
	KM(k1, k2) = ...
	    kernelStruct.kernel(kernelStruct, sts{k1}, sts2{k2}, ksize);
    end
end

%% Normalize the kernel if it has the flag
if isfield(kernelStruct, 'KM') && isfield(kernelStruct, 'normalizeKM')  && kernelStruct.normalizeKM
    if nargin < 4 % symmetric
	sd = 1./sqrt(diag(KM));
	KM = bsxfun(@times, KM, sd);
	KM = bsxfun(@times, KM', sd)';
    else
	d1 = zeros(N, 1);
	for k1 = 1:N
	    d1(k1) = kernelStruct.kernel(kernelStruct, sts{k1}, sts{k1}, ksize);
	end
	d2 = zeros(N2, 1);
	for k2 = 1:N2
	    d2(k2) = kernelStruct.kernel(kernelStruct, ...
					sts2{k2}, sts2{k2}, ksize);
	end
	d1 = 1./sqrt(d1);
	d2 = 1./sqrt(d2);
	KM = bsxfun(@times, KM, d1);
	KM = bsxfun(@times, KM', d2)';
    end
end
