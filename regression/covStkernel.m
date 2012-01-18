function K = covStkernel(metaParam, hyp, x, z, i)
% Flexible interface for GPML library (v3.1) to use spike train kernels.
% Input
%   metaParam: structure passed along from gp via feval
%	metaParam.ks spike train kernel strucuture (see kernelFactory)
%	metaParam.usePrecomputation
%	    0: recompute everything fresh
%	    1: x did not change, only recompute what involves z
%	    2: x and z did not change
%	    note that hyperparameters may have changed
%   hyp: hyper parameter vector
%	All the hyperparemters are assumed to be strictly positive.
%	A log transformation is done to keep it that way with unconstrained
%	optimization.
%	The last of the vector is the log signal variance (sf2)
%   x: spike trains data
%   z: 'diag' or spike train data
%   i: flag to retrieve derivative
%
% Behavior:
% Case covTest(metaParam) = # of hyper-parameters (use eval to get the number)
% Case covTest(metaParam, hyp, x) = K(x, x; hyp)
% Case covTest(metaParam, hyp, x, 'diag') = K(x, x; hyp)  (diagonal only)
% Case covTest(metaParam, hyp, x, z) = K(x, z; hyp)  (cross-covariance)
% Case covTest(metaParam, hyp, x, z, i = k) derivative of K(x,z;hyp)  wrt hyp(k)

persistent p_Kxz p_Kxx p_dKxz p_dKxx p_hyp

ks = metaParam.ks;
% fprintf('covStkernel %s [%d][%d]\n', ks.name, nargin, nargout);

if nargin < 2; K = sprintf('%d', ks.nParams + 1); return; end
if nargin < 3; error('Need input!'); end
if nargin > 3 && ischar(z) && strcmp(z, 'diag'); rDiag = true;
else rDiag = false; end

if nargin < 5; returnDerivative = false; else returnDerivative = true; end

assert(length(hyp) == ks.nParams + 1, '# hyperparameter includes an additional scale at the end');
assert(iscell(x), 'x needs to be a cell array of spike trains');
nX = length(x);
hyp(1:end-1) = exp(hyp(1:end-1));
% fprintf('%f\t', hyp); fprintf('\n');
sf2 = exp(hyp(end)); % signal variance hyperparameter
if rDiag && ~returnDerivative
    % we only need the diagonal of K(x, x).
    if isfield(ks, 'isShiftInvariant') && ks.isShiftInvariant
	K = sf2 * ks.kernel(ks, [], [], hyp) * ones(nX, 1);
    else
	for k = 1:nX
	    K = sf2 * ks.kernel(ks, x{k}, x{k}, hyp);
	end
	K = diag(K);
    end
    return;
end

if nargin > 3 && ~isempty(z)
    assert(iscell(z), 'z needs to be a cell array of spike trains');
else
    z = [];
end

if isfield(metaParam, 'usePrecomputation') && metaParam.usePrecomputation ~= 0
    %fprintf('Using precomputed kernel matrix or its derivatives\n');
    if isempty(p_hyp) || ~all(p_hyp == hyp)
	%fprintf('But! The hyperparameter changed, so I cannot do this\n');
	metaParam.usePrecomputation = 0;
	p_Kxx = []; p_Kxz = []; p_dKxx = {}; p_dKxz = {};
    end
else
    metaParam.usePrecomputation = 0;
    p_Kxx = []; p_Kxz = []; p_dKxx = {}; p_dKxz = {};
end

if ~returnDerivative
    if ~isempty(z)
	if metaParam.usePrecomputation >= 2 && ~isempty(p_Kxz)
	    K = p_Kxz;
	else
	    K = sf2 * computeKernelMatrix(ks, x, hyp, z);
	    p_Kxz = K;
	end
    else
	if metaParam.usePrecomputation >= 1 && ~isempty(p_Kxx)
	    K = p_Kxx;
	else
	    K = sf2 * computeKernelMatrix(ks, x, hyp);
	    p_Kxx = K;
	end
    end
else
    % derivatives
    if i == ks.nParams + 1
	if ~isempty(z)
	    if metaParam.usePrecomputation >= 2 && ~isempty(p_dKxz) && length(p_dKxz) >= i
		K = p_dKxz{i};
	    else
		K = sf2 * computeKernelMatrix(ks, x, hyp, z);
		p_dKxz{i} = K;
	    end
	else
	    if metaParam.usePrecomputation >= 1 && ~isempty(p_dKxx) && length(p_dKxx) >= i
		K = p_dKxx{i};
	    else
		K = sf2 * computeKernelMatrix(ks, x, hyp);
		p_dKxx{i} = K;
	    end
	end
    else
	if ~isempty(z)
	    if metaParam.usePrecomputation >= 2 && ~isempty(p_dKxz) && length(p_dKxz) >= i
		K = p_dKxz{i};
	    else
		K = sf2 * computeDKernelMatrix(ks, x, hyp, z);
		K = squeeze(K(:,:,i)) * hyp(i); % chain rule
		p_dKxz{i} = K;
	    end
	else
	    if metaParam.usePrecomputation >= 1 && ~isempty(p_dKxx) && length(p_dKxx) >= i
		K = p_dKxx{i};
	    else
		K = sf2 * computeDKernelMatrix(ks, x, hyp);
		K = squeeze(K(:,:,i)) * hyp(i); % chain rule
		p_dKxx{i} = K;
	    end
	end
    end
end

p_hyp = hyp;

if rDiag
    K = diag(K);
end
