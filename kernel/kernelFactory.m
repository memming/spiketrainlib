function [kernelStruct] = kernelFactory(kernelName, T, param1, param2)
% Return a spike train kernel structure.
%
% Examples:
% nci2 = kernelFactory('nCI2', 3, 'gaussian');
%
% Copyright 2010-2012 Memming. All rights reserved.

if nargin < 4
    param2 = [];
end

if nargin < 3
    param1 = [];
end

%% preprocessing kernel name for extra parameters
switch lower(kernelName)
    case {'stratified_poisson_scaling', 'stratified poisson scaling'}
	kernelName = 'stratified';
	param2 = 'poissonscaling';
end

kernelStruct.T = T;
switch lower(kernelName)
    case {'count'}
	% Binned linear kernel
	kernelStruct.name = 'Total Count';
	kernelStruct.isPSD = false;
	kernelStruct.isHermitian = true;
	kernelStruct.isSPD = false;
	kernelStruct.isShiftInvariant = true;
        kernelStruct.kernel = @(ks, st1, st2, p) (length(st1) * length(st2));
	kernelStruct.nParams = 1;
        kernelStruct.autoParam = @(ks,sts)({1});
        kernelStruct.dkernel = @(ks,x,y,ksize)(0);
    case {'count-bias'}
	% Binned linear kernel
	kernelStruct.name = 'Total Count with bias';
	kernelStruct.isPSD = false;
	kernelStruct.isHermitian = true;
	kernelStruct.isSPD = false;
	kernelStruct.isShiftInvariant = true;
        kernelStruct.kernel = @(ks, st1, st2, p) (length(st1) * length(st2) + 1);
	kernelStruct.nParams = 0;
        kernelStruct.autoParam = @(ks,sts)({1});
    case {'stratified_delta'}
	% Binned linear kernel
	kernelStruct.name = 'Stratified Delta';
	kernelStruct.isPSD = false;
	kernelStruct.isHermitian = true;
	kernelStruct.isSPD = false;
	kernelStruct.isShiftInvariant = true;
        kernelStruct.kernel = @(ks, st1, st2, p) (length(st1) == length(st2));
	kernelStruct.nParams = 0;
        kernelStruct.autoParam = @(ks,sts)({1});
    case {'blin'}
	% Binned linear kernel
	kernelStruct.name = 'Binned Linear';
	kernelStruct.isPSD = false;
	kernelStruct.isHermitian = true;
	kernelStruct.isSPD = false;
	kernelStruct.isShiftInvariant = true;
        kernelStruct.kernel = @binnedLinear;
	kernelStruct.nParams = 1;
        kernelStruct.autoParam = @(ks,sts)({ks.T / 5; ks.T / 10; ks.T / 20});
        % kernelStruct.autoParam = @(ks,sts)(autoKernelSizeScalar(@pairwiseL1, sts));
    case {'brbf'}
	% Binned linear kernel
	kernelStruct.name = 'Binned RBF';
	kernelStruct.isPSD = false;
	kernelStruct.isHermitian = true;
	kernelStruct.isSPD = false;
	kernelStruct.isShiftInvariant = true;
        kernelStruct.kernel = @binnedRBF;
	warning('BRBF-Implementtion incomplete');
	kernelStruct.nParams = 2;
        kernelStruct.autoParam = @(ks,sts) autoKernelSize_BRBF(ks, sts);
    case {'ci', 'cross intensity', 'mci', 'memoryless ci'}
	kernelStruct.name = 'mCI';
	kernelStruct.isPSD = true;
	kernelStruct.isHermitian = true;
	kernelStruct.isSPD = false;
	kernelStruct.isShiftInvariant = true;
        kernelStruct.kernel = @mci;
        kernelStruct.dkernel = @dmci;
	kernelStruct.nParams = 1;
        kernelStruct.autoParam = @(ks,sts)(autoKernelSizeScalar(@pairwiseL1, sts));
    case {'wmci', 'windowed mci'}
	kernelStruct.name = 'wmCI';
	kernelStruct.isPSD = true;
	kernelStruct.isHermitian = true;
	kernelStruct.isSPD = false;
	kernelStruct.isShiftInvariant = true;
	kernelStruct.isContinuousShift = true;
        kernelStruct.kernel = @wmci;
        kernelStruct.dkernel = @wdmci;
	kernelStruct.nParams = 1;
	[wfh, dwfh] = windowFunctionFactory(T, param1);
	kernelStruct.windowFunction = wfh;
	kernelStruct.dwindowFunction = dwfh;
        kernelStruct.autoParam = @(ks,sts)(autoKernelSizeScalar(@pairwiseL1, sts));
    case {'nci1', 'schoenberg', 'sch'}
	kernelStruct.name = ['sch (' param1 ')'];
	kernelStruct.isPSD = true;
	kernelStruct.isHermitian = true;
	kernelStruct.isSPD = true;
	kernelStruct.isShiftInvariant = true;
        kernelStruct.kernel = @nci1;
        kernelStruct.KM = @nci1_KM;
        kernelStruct.K = scalarKernelFactory(param1);
        kernelStruct.autoParam = @(ks,sts)(autoKernelSize2(...
            @pairwiseL1, @nci1_autoInner, sts));
    case {'gnci1'}
        kernelStruct.innerKS = param2;
        kernelStruct.name = ['gnCI1 [' kernelStruct.innerKS.name ']'];
        kernelStruct.isPSD = true;
        kernelStruct.isHermitian = true;
        kernelStruct.isSPD = false;
        kernelStruct.isShiftInvariant = true;
        kernelStruct.kernel = @gnci1;
        kernelStruct.KM = @gnci1_KM;
        kernelStruct.K = scalarKernelFactory(param1);
        kernelStruct.autoParam = @(ks,sts)(autoKernelSize2(...
            @pairwiseL1, @(st1,st2,ksize)gnci1_autoInner(ks,st1,st2,ksize), sts));
        kernelStruct.nParams = 2;
    case {'nci2'}
	kernelStruct.name = 'nCI2';
	kernelStruct.isPSD = true;
	kernelStruct.isHermitian = true;
	kernelStruct.isSPD = false;
	kernelStruct.isShiftInvariant = true;
        kernelStruct.kernel = @nci2;
	kernelStruct.K = scalarKernelFactory(param1);
        kernelStruct.autoParam = @(ks,sts)(autoKernelSize2(...
	    @pairwiseL1, @nci2_autoInner, sts));
    case {'i_k_int', 'i-k-int'}
	kernelStruct.name = 'I-K-int';
	kernelStruct.isPSD = true;
	kernelStruct.isHermitian = true;
	kernelStruct.isSPD = true;
	kernelStruct.isShiftInvariant = true;
        kernelStruct.kernel = @I_K_int;
	kernelStruct.K = scalarKernelFactory(param1);
        kernelStruct.autoParam = @(ks,sts)(autoKernelSizeScalar(@(x,y)I_exp_int_inner(x,y,T), sts));
    case {'i_exp_int', 'i-exp-int'}
	kernelStruct.name = 'I-exp-int';
	kernelStruct.isPSD = true;
	kernelStruct.isHermitian = true;
	kernelStruct.isSPD = true;
	kernelStruct.isShiftInvariant = true;
        kernelStruct.kernel = @I_exp_int;
        kernelStruct.dkernel = @dI_exp_int;
	kernelStruct.nParams = 1;
        kernelStruct.autoParam = @(ks,sts)(autoKernelSizeScalar(@(x,y)I_exp_int_inner(x,y,T), sts));
    case {'wi_exp_int', 'wi-exp-int'}
	kernelStruct.name = 'wI-exp-int';
	kernelStruct.isPSD = true;
	kernelStruct.isHermitian = true;
	kernelStruct.isSPD = true;
	kernelStruct.isShiftInvariant = true;
        kernelStruct.kernel = @wI_exp_int;
        kernelStruct.dkernel = @dwI_exp_int;
	[wfh, dwfh] = windowFunctionFactory(T, param1);
	kernelStruct.windowFunction = wfh;
	kernelStruct.dwindowFunction = dwfh;
	kernelStruct.nParams = 1;
        kernelStruct.autoParam = @(ks,sts)(autoKernelSizeScalar(@(x,y)wI_exp_int_inner(x,y,T,ks.windowFunction), sts));
    case {'i_int_exp', 'i-int-exp'}
	kernelStruct.name = 'I-int-exp';
	kernelStruct.isPSD = true;
	kernelStruct.isHermitian = true;
	kernelStruct.isSPD = true;
	kernelStruct.isShiftInvariant = true;
        kernelStruct.kernel = @I_int_exp;
        kernelStruct.autoParam = @(ks,sts)(autoKernelSizeScalar(@(x,y)I_int_exp_inner(x,y,T), sts));
    case {'reek', 'nfisher', 'reef', 'reefkernel'} 
	% Nicholas (Neko) Fisher, A. Banerjee NIPS 2010
	kernelStruct.name = 'REEF kernel';
	kernelStruct.desc = '(Fisher and Banerjee, 2010)';
	kernelStruct.isPSD = true;
	kernelStruct.isHermitian = true;
	kernelStruct.isSPD = NaN;
	kernelStruct.isShiftInvariant = false;
        kernelStruct.kernel = @reek;
        %kernelStruct.autoParam = @(sts)(max(stsfun(@max,sts)));
        kernelStruct.autoParam = @(ks,sts)({1});
    case {'schrauwen'}
	% Benjamin Schrauwen's kernels, mCI is a special case with Laplacian
	kernelStruct.name = 'Schrauwen';
	kernelStruct.desc = '(Schrauwen and Campenhout, 2007)';
	kernelStruct.isPSD = true;
	kernelStruct.isHermitian = true;
	kernelStruct.isSPD = NaN; % not true..I think
	kernelStruct.isShiftInvariant = true;
        kernelStruct.kernel = @schrauwen;
	kernelStruct.K = scalarKernelFactory(param1);
        kernelStruct.autoParam = @(ks,sts)(autoKernelSizeScalar(@pairwiseL1, sts));
    case {'stratified'}
	% Stratified kernel
	kernelStruct.isPSD = true;
	kernelStruct.isHermitian = true;
	kernelStruct.isSPD = true; % when K is
	kernelStruct.isShiftInvariant = true;
        kernelStruct.kernel = @stratified;
	kernelStruct.K = scalarKernelFactory(param1);
	if strcmp(lower(param2), 'poissonscaling')
	    kernelStruct.name = 'Stratified (PS)';
	    kernelStruct.isPoissonScaling = true;
	    kernelStruct.scalingFactor = [0.634401844206212 0.787921638254979 0.854916422727872 0.892257481024682 0.916198544928009 0.931426092471907 0.942631790652277 0.952587333162107 0.958578082280984 0.964914744588260 0.969725608428331 0.973847376439812 0.976107306928358 0.980130879195747 0.981848557078916 0.983951004809557 0.986104485016542 0.988075672059994 0.988959169321910 0.990977825008080 0.991542854220363 0.993225659560109 0.994530755867158 0.995324788662448 0.996065963841616 0.996811773289903 0.997760611096988 0.999153476479780 0.999362983762747 1.000000000000000];
	else
	    kernelStruct.name = 'Stratified (constant)';
	    kernelStruct.isPoissonScaling = false;
	end
        kernelStruct.autoParam = @(ks,sts)(autoKernelSizeScalar_stratified(ks, @stratifiedL2, sts));
    case {'stratified_min'}
	% Stratified kernel
	kernelStruct.name = 'Stratified min';
	kernelStruct.isPSD = true;
	kernelStruct.isHermitian = true;
	kernelStruct.isSPD = true; % when K is
	kernelStruct.isShiftInvariant = false;
        kernelStruct.kernel = @stratified_min;
	kernelStruct.K = scalarKernelFactory(param1);
        kernelStruct.autoParam = @(ks,sts)({1});
    case {'spikernel'} % Spikernel
	%   Spikernels: Embedding Spiking Neurons in Inner-Product Spaces
	%   Lavi Shpigelman, Yoram Singer, Rony Paz and Eilon Vaadia
	%   Advances in Neural Information Processing Systems (NIPS) 15
	%   MIT Press, Cambridge, MA, 2003.
	
	if ~exist('Fspike.m')
	    kernelStruct = [];
	    return
	end
	kernelStruct.name = 'spikernel';
	kernelStruct.desc = '(Shpigelman et al. NC 2003)';
	kernelStruct.isPSD = true;
	kernelStruct.isHermitian = true;
	kernelStruct.isSPD = false; % not known
	kernelStruct.isShiftInvariant = true;
        kernelStruct.kernel = @spikernel;
        kernelStruct.autoParam = @autoKernelSize_spikernel;
        kernelStruct.paramSpec = {'bin size', 'N', 'lam', 'mu', 'p'};
    case {'alignment'} % Alignment kernel (TODO)
	warning('Alignment kernel is not implemented');
	kernelStruct.kernel = @(x,y,z)([]);
        kernelStruct.autoParam = @(ks,sts)({});
	kernelStruct.desc = '(Eichhorn et al. NIPS 2003)';
    case {'metric-pseudo', 'seth2011a'}
	% Sohan Seth's SPD kernel derived from arbitrary spike train metric
	kernelStruct.name = 'metric based non-PD kernel';
	kernelStruct.isPSD = false;
	kernelStruct.isHermitian = true;
	kernelStruct.isSPD = false;
        kernelStruct.kernel = @metricPseudoKernel;
	kernelStruct.stmetric = metricFactory(param1);
	kernelStruct.cmf = cmfFactory(param2);
        kernelStruct.autoParam = @(ks,sts)(autoKernelSizeScalar(@kernelStruct.stmetric, sts));
    case {'spike-align-kernel', 'seth2011b'}
	kernelStruct.name = 'Spike alignment kernel';
	kernelStruct.isPSD = 0;
	kernelStruct.isHermitian = true;
	kernelStruct.isSPD = 0;
        kernelStruct.kernel = @spikeAlignKernel;
        kernelStruct.autoParam = @(ks,sts)(autoKernelSizeScalar(@pairwiseL1, sts));
    case { 'normalized-spike-align-kernel', 'normalized-seth2011b', 'n-seth2011b'}
	kernelStruct = getNormalizedKernel('seth2011b', T, param1, param2);
    	kernelStruct.name = 'Normalized spike alignment kernel';
    case {'spike-align-kernel-isi', 'spike-align-kernel-2', 'seth2011c'}
	kernelStruct.name = 'Spike alignment kernel ISI';
	kernelStruct.isPSD = 0;
	kernelStruct.isHermitian = true;
	kernelStruct.isSPD = 0;
        kernelStruct.kernel = @spikeAlignKernel2;
        kernelStruct.nParams = 1;
	% TODO pairwiseL1 could be all zeros => Inf kernel size problem
        kernelStruct.autoParam = @(ks,sts)(autoKernelSizeScalar(@(x,y)(1./(eps+pairwiseL1(x,y))), sts));
    case { 'normalized-spike-align-kernel-2', 'normalized-spike-align-kernel-isi', 'normalized-seth2011c', 'n-seth2011c'}
	kernelStruct = getNormalizedKernel('seth2011c', T, param1, param2);
    	kernelStruct.name = 'Normalized spike alignment kernel ISI';
    case {'seth2011c'}
	kernelStruct.name = 'TEMP';
	kernelStruct.desc = 'TEMP';
	kernelStruct.isPSD = true;
	kernelStruct.isHermitian = true;
	kernelStruct.isSPD = false;
	kernelStruct.isShiftInvariant = false;
        kernelStruct.kernel = @inhomogeneousPoissonLogitNormalKernel;
        kernelStruct.autoParam = @(ks,sts)({1});
    otherwise
        error('Unknown kernel name [%s]', kernelName);
end
kernelStruct.isMercer = and(kernelStruct.isPSD, kernelStruct.isHermitian);

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [kernelStruct] = getNormalizedKernel(kernelName, T, param1, param2)
% normalize a kernel such that it has 1s in the diagonal
% k(x,y) = k'(x,y) / sqrt( k(x,x) * k(y,y) )
kernelStruct = kernelFactory(kernelName, T, param1, param2);
kOriginal = kernelStruct.kernel;
kernelStruct.kernel = @(ks,x,y,ksize) (kOriginal(ks,x,y,ksize) / sqrt(kOriginal(ks,x,x,ksize) * kOriginal(ks,y,y,ksize)));
kernelStruct.normalizeKM = true;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function K = scalarKernelFactory(kern_str)

if nargin < 1 || isempty(kern_str)
    kern_str = 'gaussian';
end

switch lower(kern_str)
case 'gaussian'
    K = @(z,ksize) exp(-(z.^2) ./ (2*(ksize^2)));
case 'laplacian'
    K = @(z,ksize) exp(-abs(z) ./ ksize);
case 'triangular'
    K = @(z,ksize) ((abs(z) < ksize) .* (1 - abs(z)/ksize));
case 'rectwin'
    K = @(z,ksize) (abs(z) <= ksize);
otherwise
    error('Unknown kernel! Try one: ''laplacian'', ''gaussian'', ''triangular'' and ''rectwin''');
end
end % end sub function scalarKernelFactory

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function cmf = cmfFactory(cmfName)
% Ref: Miller and Samko, Completely monotonic functions, 
%      Integr. Transf. and Spec. Funct. 2001, v. 12, n. 4, 389--402

if nargin < 1 || isempty(cmfName)
    cmfName = 'laplacian';
end

switch lower(cmfName)
case 'laplacian'
    cmf = @(x) exp(-x);
case 'cauchylike'
    lambda = 1;
    nu = 2;
    cmf = @(x) 1./(lambda + x).^nu; % lambda >= 0 nu >= 0
case 'lninv'
    b = 1;
    cmf = @(x) log(b + 1/x); % b >= 1
case 'invexp'
    cmf = @(x) exp(1 ./ x); % (Miller & Samko 2001, eq 1.3)
case 'lnoverx'
    cmf = @(x) ln(1 + x) ./ x; % (Miller & Samko 2001, eq 1.3)
case 'explike'
    a = 1;
    cmf = @(x) (1 + x).^(a./x); % a >= 0 (Miller & Samko 2001, eq 1.18)
end
end % sub function cmfFactory

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function v = nci1(ks, st1, st2, ksize)
    m1 = mci([], st1, st1, ksize(1));
    m2 = mci([], st1, st2, ksize(1));
    m3 = mci([], st2, st2, ksize(1));
    D2 = m1 + m3 - 2*m2;
    v = ks.K(sqrt(D2), ksize(2));
end

function D = nci1_autoInner(st1, st2, ksize1)
    m1 = mci([], st1, st1, ksize1);
    m2 = mci([], st1, st2, ksize1);
    m3 = mci([], st2, st2, ksize1);
    D = sqrt(abs(m1 + m3 - 2*m2));
end

function mKM = nci1_KM(ks, sts, ksize)
    mciKS = kernelFactory('mci', ks.T);
    mKM = computeKernelMatrix(mciKS, sts, ksize(1));
    g = diag(mKM);
    o = ones(1, size(mKM,1));
    mKM = g(:) * o + o' * g(:)' - 2 * mKM; % square distance
    mKM = ks.K(sqrt(mKM), ksize(2)); % same variable to save memory
end

%% Generalized nCI1 with arbitrary kernel
function v = gnci1(ks, st1, st2, ksize)
    m1 = ks.innerKS.kernel(ks.innerKS, st1, st1, ksize(1));
    m2 = ks.innerKS.kernel(ks.innerKS, st1, st2, ksize(1));
    m3 = ks.innerKS.kernel(ks.innerKS, st2, st2, ksize(1));
    D2 = m1 + m3 - 2*m2;
    v = ks.K(sqrt(D2), ksize(2));
end

function D = gnci1_autoInner(ks, st1, st2, ksize1)
    m1 = ks.innerKS.kernel(ks.innerKS, st1, st1, ksize1);
    m2 = ks.innerKS.kernel(ks.innerKS, st1, st2, ksize1);
    m3 = ks.innerKS.kernel(ks.innerKS, st2, st2, ksize1);
    D = sqrt(abs(m1 + m3 - 2*m2));
end

function mKM = gnci1_KM(ks, sts, ksize)
    mKM = computeKernelMatrix(ks.innerKS, sts, ksize(1));
    g = diag(mKM);
    o = ones(1, size(mKM,1));
    mKM = g(:) * o + o' * g(:)' - 2 * mKM; % square distance
    mKM = ks.K(sqrt(mKM), ksize(2)); % same variable to save memory
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function v = nci2(ks, st1, st2, ksizes)
% Evaluates the nonlinear Cross-Intensity kernel (definition 2 in Paiva 2009).
% \int K_ksize(\lambda1(t) - \lambda2(t)) dt
% where lambda is the smoothed spike train with a rectangular function of pwidth
pwidth = ksizes(1);
ksize = ksizes(2);
v = 0;
if isempty(st1) && isempty(st2)
    v = ks.K(0, ksize);
    return;
end

[val, times] = nci2_inner(st1, st2, pwidth);

times(times < 0) = 0;
times(times > ks.T) = ks.T;
dtimes = diff([times ks.T]);
nzidx = (dtimes ~= 0);
v = sum(dtimes(nzidx) .* ks.K(val(nzidx), ksize)) / ks.T; % integrate over time
end % sub function nci2

function [val, times] = nci2_inner(st1, st2, pwidth)
st1 = st1(:)'; st2 = st2(:)';
L1 = length(st1); L2 = length(st2);
% times: when the difference between intensity functions changes
[times idx] = sort([0, st1-pwidth/2, st1+pwidth/2, st2-pwidth/2, st2+pwidth/2]);
% if the difference should increase or decrease
incr = [0, ones(1,L1), -ones(1,L1), -ones(1,L2), ones(1,L2)];
incr = incr(idx);
val = cumsum(incr) / pwidth;
end % sub function nci2_inner

function val = nci2_autoInner(st1, st2, pwidth)
val = nci2_inner(st1, st2, pwidth);
val = val(val~=0);
val = abs(val);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function v = mci(ks, st1, st2, ksize)
% memoryless cross intensity kernel (mCI)
% Input:
%   st1: (N1x1) sorted spike times
%   st2: (N2x1) sorted spike times
%   ksize: kernel size
% Output:
%   v: sum_i sum_j exp(-|st1{i] - st{j}|/sigma)

N1 = length(st1);
N2 = length(st2);

v = 0;
if N1 == 0 || N2 == 0; return; end

if ~(isreal(ksize(1)) && ksize(1) > 0)
    error('Kernel size must be non-negative real');
end

if N1 <= 4 || N2 <= 4
    inside = pairwiseL1(st1, st2);
    v = sum(exp(-inside(:)/ksize(1)));
else
    % O(N log(N)) implementation
    X = st1(:) / ksize(1); Y = st2(:) / ksize(1);
    X_sort = sort(X); Y_sort = sort(Y);

    X_pos_exp = exp(X_sort); X_neg_exp = exp(-X_sort);
    Y_pos_exp = exp(Y_sort); Y_neg_exp = exp(-Y_sort);

    Y_pos_cum_sum_exp  = cumsum(Y_pos_exp);
    Y_neg_cum_sum_exp  = flipud(cumsum(flipud(Y_neg_exp)));

    Y_sort = [Y_sort; Inf]; yidx = 0;

    for xidx = 1:length(X)
	while Y_sort(yidx+1) <= X_sort(xidx); yidx = yidx + 1; end

	if yidx == 0
	    v = v + Y_neg_cum_sum_exp(1) * X_pos_exp(xidx);
	elseif yidx == length(Y)
	    v = v + Y_pos_cum_sum_exp(end) * X_neg_exp(xidx);
	else
	    v = v + Y_pos_cum_sum_exp(yidx) * X_neg_exp(xidx) ...
		+ Y_neg_cum_sum_exp(yidx+1) * X_pos_exp(xidx);
	end
    end
end
end % end mci

function v = dmci(ks, st1, st2, ksize)
    inside = pairwiseL1(st1, st2);
    v = sum(inside(:) .* exp(-inside(:)/ksize(1))) / ksize(1)^2;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function v = wmci(ks, st1, st2, ksize)
% windowed memoryless cross intensity kernel (mCI)
% Note the higher time complexity.
%
% Input:
%   st1: (N1x1) sorted spike times
%   st2: (N2x1) sorted spike times
%   ksize: kernel size
% Output:
%   v: sum_i sum_j exp(-|st1{i] - st{j}|/sigma) * f_{window}(st{i}) * f_{window}(st{j})

N1 = length(st1);
N2 = length(st2);

v = 0;
if N1 == 0 || N2 == 0; return; end

w1 = ks.windowFunction(st1);
w2 = ks.windowFunction(st2);

if ~(isreal(ksize(1)) && ksize(1) > 0)
    error('Kernel size must be non-negative real');
end

inside = pairwiseL1(st1, st2);
v = exp(-inside / ksize(1)) .* (w1(:) * w2(:)');
v = sum(v(:));

end % end wmci

function v = dwmci(ks, st1, st2, ksize)
    inside = pairwiseL1(st1, st2);
    w1 = ks.windowFunction(st1);
    dw1 = ks.dwindowFunction(st1);
    w2 = ks.windowFunction(st2);
    dw2 = ks.dwindowFunction(st1);
    v = exp(-inside / ksize(1)) .* (w1(:) * w2(:)');
    v = v .* (inside / ksize(1)^2 + repmat(dw1(:), 1, length(st2)) + repmat(dw2(:)', length(st2), 1));
    v = sum(v);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function d = pairwiseL1(x, y)

d = bsxfun(@minus, repmat(x(:), 1, length(y)), y(:)');
d = abs(d);

% d = pdist2(x(:), y(:), 'minkowski', 1);
end % pairwiseL1

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [z, zz] = I_common(st1, st2)
    l1 = length(st1) + 1; % empty spike train
    [z, cls] = sort([st1(:); 0; st2(:); 0]);
    zz = ones(size(z));
    zz(cls > l1) = -1;
    zz = cumsum(zz);
    % alternatively
    % zz2 = [ones(numel(st1),1); 1; -ones(numel(st2),1); -1];
    % zz2 = zz2(cls); zz2 = cumsum(zz2);

    % alternatively
    % z1 = cumsum(cls <= l1);
    % z2 = cumsum(cls > l1);
    % zz = (z1 - z2);
    zz = zz.^2;
end

function [kk] = I_identity(st1, st2, T) % This is not a kernel (not used)
    [z, zz] = I_common(st1, st2);
    kk = sum(diff([z; T]) .* zz);
end

function [zzz] = I_exp_int_inner(st1, st2, T)
    [z, zz] = I_common(st1, st2);
    zzz = sum(diff([z; T]) .* zz); % note that zz is a squared quantity
end

function [kk] = I_exp_int(ks, st1, st2, ksize)
    [zzz] = I_exp_int_inner(st1, st2, ks.T);
    kk = exp(-zzz/ksize(1));
end

function [kk] = dI_exp_int(ks, st1, st2, ksize)
    [zzz] = I_exp_int_inner(st1, st2, ks.T);
    kk = exp(-zzz/ksize(1)) * zzz / ksize(1)^2;
end

function [kk] = I_K_int(ks, st1, st2, ksize)
    [z, zz] = I_common(st1, st2);
    zzz = sum(diff([z; ks.T]) .* zz);
    kk = ks.K(sqrt(zzz),ksize);
end

function [zz] = I_int_exp_inner(st1, st2, T)
    [z, zz] = I_common(st1, st2);
    zz = zz(zz ~= 0);
end

function [kk] = I_int_exp(ks, st1, st2, ksize)
    [z, zz] = I_common(st1, st2);
    zz = exp(-zz/ksize);
    kk = sum(diff([z; ks.T]) .* zz) / ks.T;
end

%% Windowed I kernels
function [z, zz] = wI_common(st1, st2, windowFunction)
    l1 = length(st1) + 1; % empty spike train
    [z, cls] = sort([st1(:); 0; st2(:); 0]);
    zz = [windowFunction(st1(:)); 0; -windowFunction(st2(:)); 0];
    zz = zz(cls);
    zz = cumsum(zz);
    zz = zz.^2;
end

function [zzz] = wI_exp_int_inner(st1, st2, T, windowFunction)
    [z, zz] = wI_common(st1, st2, windowFunction);
    zzz = sum(diff([z; T]) .* zz); % note that zz is a squared quantity
end

function [kk] = wI_exp_int(ks, st1, st2, ksize)
    [zzz] = wI_exp_int_inner(st1, st2, ks.T, ks.windowFunction);
    kk = exp(-zzz/ksize(1));
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function v = reek(ks, st1, st2, ksize)
    v = 0;
    st1 = ks.T - st1;
    st2 = ks.T - st2;
    for k1 = 1:length(st1)
	% REEF kernel is only defined for t ~= 0
	if st1(k1) == 0; continue; end
	for k2 = 1:length(st2)
	    if st2(k2) == 0; continue; end
	    v = v + (st1(k1) * st2(k2)) / (st1(k1) + st2(k2))^2;
	end
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function candidates = autoKernelSize_BRBF(ks, spikeTrains)
    nSpikeTrains = numel(spikeTrains);
    M = min(200, nSpikeTrains);
    candidates = cell(9,1);

    binSizeList = [ks.T/10; ks.T/20; ks.T/40]; % For 3 different bin sizes
    for k = 1:3
	k1 = randperm(nSpikeTrains); k1 = k1(1:M);
	k2 = randperm(nSpikeTrains); k2 = k2(1:M);
	% x1, x2: binlen x #sptrains
	x1 = binSpikeTrains(spikeTrains(k1), ks.T, binSizeList(k));
	x2 = binSpikeTrains(spikeTrains(k2), ks.T, binSizeList(k));
	d = sqrt(sum((x1 - x2).^2));
	d = d(d~=0);
	cand2 = quantile(d, [0.1 0.5 0.9]); % Strictly positive
	cand2 = [ cand2(1)/2 cand2 cand2(end) * 2 ];
	%cand2 = median(d) * 2.^[-2:2]'; % Median
	for kk = 1:length(cand2)
	    candidates{kk+(k-1)*3} = [binSizeList(k); cand2(kk)];
	end
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function v = binnedRBF(ks, st1, st2, ksize)
    x1 = binSpikeTrains({st1}, ks.T, ksize(1));
    x2 = binSpikeTrains({st2}, ks.T, ksize(1));
    v = exp(-sum((x1 - x2).^2)/2/ksize(2)^2);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function v = binnedLinear(ks, st1, st2, ksize)
    x1 = binSpikeTrains({st1}, ks.T, ksize);
    x2 = binSpikeTrains({st2}, ks.T, ksize);
    v = x1' * x2 + 1; % 1 is for bias
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function v = schrauwen(ks, st1, st2, ksize)
    if isempty(st1) || isempty(st2)
	v = 0;
	return;
    end
    inside = pairwiseL1(st1, st2);
    v = sum(ks.K(inside(:),ksize));
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function d = stratifiedL2(x, y)

if length(x) ~= length(y)
    d = [];
    return;
end

if isempty(x)
    d = 0;
    return;
end

d = sqrt(sum((x(:) - y(:)).^2));

end % stratifiedL2

function v = stratified(ks, st1, st2, ksize)

N1 = length(st1); N2 = length(st2);

if N1 ~= N2
    v = 0;
    return;
end

if N1 == 0
    v = ks.K(0, ksize);
    return;
end

if ks.isPoissonScaling && N1 <= length(ks.scalingFactor)
    ksize = ksize * ks.scalingFactor(N1);
end

v = ks.K(sqrt(sum((st1 - st2).^2)), ksize);

end % stratified

function candidates = autoKernelSizeScalar_stratified(ks, scaledQuantityHandle, spikeTrains)
% For stratified kernel, subsample smart such that dimensions match
nSpikeTrains = numel(spikeTrains);
NN = cellfun('length', spikeTrains);

% Shall we subsample or use all possible pairwise values?
% If we use all possible pairwise values, that would not match the natural
% statistic of random pairs. These two options are not compatible.
% % Compute the max possible non-zero/non-infinite distance pairs
% % sum(nchoosek(NN,2))

M = min(200, nSpikeTrains);
k1 = randperm(nSpikeTrains);
k1 = k1(1:M);

%% Computing the quantity for the subsampled pairs
d = [];
for k = 1:M
    N1 = length(spikeTrains{k1(k)});
    idx = find(NN == N1); % choose one with same length
    if length(idx) == 1
	continue;
    end
    k2 = idx(ceil(rand * length(idx)));
    dd = scaledQuantityHandle(spikeTrains{k1(k)}, spikeTrains{k2});
    if ks.isPoissonScaling && N1 ~= 0 && N1 <= length(ks.scalingFactor)
	dd = dd / ks.scalingFactor(N1);
    end
    d = [d; dd(:)];
end

%% Choose the quantiles
if all(d == 0)
    warning('All pairwise scaled quantity is zero!');
    candidates = {Inf}; % return something large, so that it becomes count distribution comparison
    return;
else
    d = d(d~=0);
    %candidates = median(d) * 2.^[-2:2]'; % Median
    candidates = quantile(d, [0.1; 0.5; 0.9]); % Strictly positive
    candidates = [ candidates(1)/2; candidates; candidates(end) * 2 ];
    % Add twice the maximum value as well. This allows to approximate the count
    % distribution based kernel.
    %candidates(end+1) = max(d) * 2;
    % candidates = unique(candidates(:));
    %candidates = candidates(:);
end

candidates = mat2cell(candidates, ones(size(candidates,1),1), size(candidates,2));
end % autoKernelSizeScalar_stratified
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function v = stratified_min(ks, st1, st2, ksize)

if length(st1) ~= length(st2)
    v = 0;
    return;
end

if isempty(st1)
    v = 1;
    return;
end

% prod_{i=1}^d (T-max(w_1^i,w_2^i))
v = prod(ks.T-max(st1, st2));

end % stratified_min

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function v = spikernel(ks, st1, st2, param)
    binsize = param(1);
    N = param(2);
    lam = param(3);
    mu = param(4);
    p = param(5);
    binEdges = 0:binsize:ks.T;
    if isempty(st1)
	seqS = zeros(size(binEdges));
    else
	seqS = histc(st1, binEdges);
    end
    if isempty(st2)
	seqT = zeros(size(binEdges));
    else
	seqT = histc(st2, binEdges);
    end
    v = Fspikernel(seqS(:), seqT(:), N, lam, mu, p);
end

function candidates = autoKernelSize_spikernel(ks, sts)
    binSize = ks.T / 20;
    candidates = {
	[binSize,  3, 0.5, 0.7, 1]; ...
	[binSize,  5, 0.5, 0.7, 1]; ...
	[binSize, 10, 0.5, 0.7, 1]; ...
	[binSize,  3, 0.9, 0.7, 1]; ...
	[binSize,  5, 0.9, 0.7, 1]; ...
	[binSize, 10, 0.9, 0.7, 1]; ...
	[binSize,  3, 0.5, 0.9, 1]; ...
	[binSize,  5, 0.5, 0.9, 1]; ...
	[binSize, 10, 0.5, 0.9, 1]; ...
	[binSize,  3, 0.9, 0.9, 1]; ...
	[binSize,  5, 0.9, 0.9, 1]; ...
	[binSize, 10, 0.9, 0.9, 1]; ...
	}; % TODO
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function v = metricPseudoKernel(ks, st1, st2, scaleParam)
    v = ks.stmetric(st1, st2);
    v = ks.cmf(v / scaleParam);
end
