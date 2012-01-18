disp('$Id: testKernelFactory.m 92 2011-09-07 22:36:23Z memming $');

sts = cell(5, 1);
sts{1} = [];
sts{2} = [0];
sts{3} = [1];
sts{4} = [0 1];
sts{5} = [0 1 2];

T = 3;
tol = 10 * eps;

assertRange = @(x,y,msg) assert(abs(x-y) < tol, msg);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
ks = kernelFactory('brbf', T);
ksizeList = ks.autoParam(ks, sts);
for cidx = 1:numel(ksizeList)
    K = computeKernelMatrix(ks, sts, ksizeList{cidx})
end
K2 = computeKernelMatrix(ks, sts, ksizeList{cidx});
assert(abs(max(K(:) - K2(:))) == 0, 'special kernel matrix computation');

assert(ks.kernel(ks,[], [], [2 2]) == 1, 'translation invariant with diag 1');
t = rand + 0.1;
assert(ks.kernel(ks,t,t, [2 2]) == 1, 'translation invariant with diag 1');
assert(ks.kernel(ks,[t, t*2],[t t*2], [2 2]) == 1, 'translation invariant with diag 1');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
ks = kernelFactory('seth2011c', T);
ksizeList = ks.autoParam(ks, sts)
assertRange(ks.kernel(ks,[], [], 1), 1 , 'One if empty');
assertRange(ks.kernel(ks,[], [0 1 2 3], 1), 0.049787068367864, 'exp(-len) if one empty');
assertRange(ks.kernel(ks,0,2,1), mean([exp(-2) exp(-2) exp(-2)]), 'easy case1');
assertRange(ks.kernel(ks,[0 1],[0 1 2],1), 0.191161456280565, 'emp1');
assertRange(ks.kernel(ks,[0],[0 1 2],1), 0.135335283236613, 'emp2');
assertRange(ks.kernel(ks,[0 1.5],[0 1 2],1), 0.115945284189479, 'emp3');
K = computeKernelMatrix(ks, sts, ksizeList{1})
v = eig(K);
assert(all(v > -eps*5), 'Strictly positive definite, NOT?');

%% Normalized version
ks = kernelFactory('n-seth2011c', T);
ksizeList = ks.autoParam(ks, sts)
assertRange(ks.kernel(ks,[], [], 1), 1 , 'identical');
assertRange(ks.kernel(ks,[1], [1], 1), 1 , 'identical');
assertRange(ks.kernel(ks,[0 1], [0 1], 1), 1 , 'identical');
assertRange(ks.kernel(ks,0,2,1), 0.230228725493028, 'easy case1');
K = computeKernelMatrix(ks, sts, ksizeList{1})
v = eig(K);
assert(all(v > -eps*5), 'Strictly positive definite, NOT?');
assertRange(K(1,2), ks.kernel(ks,sts{1},sts{2}, ksizeList{1}), 'KM vs raw 1');
assertRange(K(4,5), ks.kernel(ks,sts{4},sts{5}, ksizeList{1}), 'KM vs raw 2');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
ks = kernelFactory('seth2011b', T);
ksizeList = ks.autoParam(ks, sts)
assertRange(ks.kernel(ks,[], [], 1), 1 , 'One if empty');
assertRange(ks.kernel(ks,[], [0 1 2 3], 1), exp(-4), 'exp(-len) if one empty');
assertRange(ks.kernel(ks,[0 1 2 3], [], 1), exp(-4), 'exp(-len) if one empty');
assertRange(ks.kernel(ks,0,2,1), mean([exp(-2) exp(-2) exp(-2)]), 'easy case1');
assertRange(ks.kernel(ks,2,0,1), mean([exp(-2) exp(-2) exp(-2)]), 'easy case1');
assertRange(ks.kernel(ks,2,1,1), mean([exp(-2) exp(-1) exp(-2)]), 'easy case2');
assertRange(ks.kernel(ks,0,1,1), mean([exp(-2) exp(-1) exp(-2)]), 'easy case2');
assertRange(ks.kernel(ks,[0 1],[0 1 2],1), 0.037873320841298, 'emp1');
assertRange(ks.kernel(ks,[0],[0 1 2],1), 0.026146525288189, 'emp2');
assertRange(ks.kernel(ks,[0 1.5],[0 1 2],1), 0.039834879992530, 'emp3');
K = computeKernelMatrix(ks, sts, ksizeList{1})
v = eig(K);
assert(all(v > 0), 'Strictly positive definite, NOT?');


%% Normalized version
ks = kernelFactory('normalized-seth2011b', T);
ksizeList = ks.autoParam(ks, sts)
assertRange(ks.kernel(ks,[], [], 1), 1 , 'identical');
assertRange(ks.kernel(ks,[1], [1], 1), 1 , 'identical');
assertRange(ks.kernel(ks,[0 1], [0 1], 1), 1 , 'identical');
K = computeKernelMatrix(ks, sts, ksizeList{1})
v = eig(K);
assert(all(v > 0), 'Strictly positive definite, NOT?');
assertRange(K(1,2), ks.kernel(ks,sts{1},sts{2}, ksizeList{1}), 'KM vs raw 1');
assertRange(K(4,5), ks.kernel(ks,sts{4},sts{5}, ksizeList{1}), 'KM vs raw 2');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
ks = kernelFactory('mci', T);
assert(abs(ks.kernel(ks,[], 0, 1)) < tol, 'Zero if empty');
assertRange(ks.kernel(ks,[], [], 1), 0 , 'Zero if empty');
assert(abs(ks.kernel(ks,0, 0, 1) - 1) < tol, 'identical single spike');
t = rand + 0.1;
assert(abs(ks.kernel(ks,t, t, 1) - 1) < tol, 'identical single spike');
assert(abs(ks.kernel(ks,0:5, 0:5, 1) - 11.146937537367616) < tol, 'identical single spike');
ksizeList = ks.autoParam(ks, sts)
K = computeKernelMatrix(ks, sts, ksizeList{1})
K2 = computeKernelMatrix(ks, sts, ksizeList{1}, sts)
assert(abs(max(K(:) - K2(:))) == 0, 'two argument compute kernel matrix');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
ks = kernelFactory('wmci', T, '1');
assert(abs(ks.kernel(ks,[], 0, 1)) < tol, 'Zero if empty');
assertRange(ks.kernel(ks,[], [], 1), 0 , 'Zero if empty');
assert(abs(ks.kernel(ks,0, 0, 1) - 1) < tol, 'identical single spike');
t = rand + 0.1;
assert(abs(ks.kernel(ks,t, t, 1) - 1) < tol, 'identical single spike');
assert(abs(ks.kernel(ks,0:5, 0:5, 1) - 11.146937537367616) < tol, 'identical single spike');

ks = kernelFactory('wmci', 4, 'sin');
assertRange(ks.kernel(ks,0:4, [1 1.4 3.9], 1), ks.kernel(ks,1:3, [1 1.4 3.9], 1), 'spikes at the boundary are ignored');
% Check continuity
assertRange(ks.kernel(ks,[1 3], [1 1.4 3.9], 1), ...
	    ks.kernel(ks,[eps 1 3], [1 1.4 3.9], 1), 'inserting before');
assertRange(ks.kernel(ks,[1 3], [1 1.4 3.9], 1), ...
	    ks.kernel(ks,[1 3 (4-eps)], [1 1.4 3.9], 1), 'inserting after');
assertRange(ks.kernel(ks,[1 3], [1 1.4 3.9], 1), ...
	    ks.kernel(ks,[1 3], [eps 1 1.4 3.9], 1), 'inserting before 2');
assertRange(ks.kernel(ks,[1 3], [1 1.4 3.9], 1), ...
	    ks.kernel(ks,[1 3], [1 1.4 3.9 (4-eps)], 1), 'inserting after 2');

ksizeList = ks.autoParam(ks, sts)
K = computeKernelMatrix(ks, sts, ksizeList{1})
K2 = computeKernelMatrix(ks, sts, ksizeList{1}, sts)
assert(abs(max(K(:) - K2(:))) == 0, 'two argument compute kernel matrix');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
ks = kernelFactory('reek', T); % REEF kernel (Fisher & Banerjee 2011)
assert(ks.kernel(ks,[], [], 0) == 0, 'zero if empty');
assert(ks.kernel(ks,[], 0, 0) == 0, 'zero if any empty');
assert(ks.kernel(ks,0, [], 0) == 0, 'zero if any empty');
assert(ks.kernel(ks,T, T, 0) == 0, 'zero if zero only');
assert(ks.kernel(ks,T, [1 2 3], 0) == 0, 'zero if zero only');
assert(ks.kernel(ks,[1 2 3], T, 0) == 0, 'zero if zero only');
assert(abs(ks.kernel(ks,[1], [1], 0) - 1/4) < tol, 'one times one over 4');
t = rand + 0.1;
assert(abs(ks.kernel(ks,t, t, 2) - 1/4) < tol, 'one times one over 4');
assert(abs(ks.kernel(ks,[1 2], [2], 3) - 1/4 - 2/9) < tol, 'one times one over 4');

ksizeList = ks.autoParam(ks, sts)
K = computeKernelMatrix(ks, sts, ksizeList{1})

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
ks = kernelFactory('nCI1', T, 'Gaussian');
ksizeList = ks.autoParam(ks, sts);
for cidx = 1:numel(ksizeList)
    K = computeKernelMatrix(ks, sts, ksizeList{cidx})
end
rmfield(ks, 'KM');
K2 = computeKernelMatrix(ks, sts, ksizeList{cidx});
assert(abs(max(K(:) - K2(:))) == 0, 'special kernel matrix computation');

assert(ks.kernel(ks,[], [], [2 2]) == 1, 'translation invariant with diag 1');
t = rand + 0.1;
assert(ks.kernel(ks,t,t, [2 2]) == 1, 'translation invariant with diag 1');
assert(ks.kernel(ks,[t, t*2],[t t*2], [2 2]) == 1, 'translation invariant with diag 1');
assert(ks.kernel(ks,0,1,[1 1]) == exp(-(2-2*exp(-1))/2), 'translation invariant with diag 1');
assert(ks.kernel(ks,0,1,[1 2]) == exp(-(2-2*exp(-1))/2/2^2), 'translation invariant with diag 1');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
ks = kernelFactory('nCI2', T, 'Gaussian');
ksizeList = ks.autoParam(ks, sts)
for cidx = 1:size(ksizeList,1)
    K = computeKernelMatrix(ks, sts, ksizeList{cidx})
end

assert(ks.kernel(ks,[], [], [1 1]) == 1, 'translation invariant w diag K(0)');
assertRange(ks.kernel(ks,1,[], [1 1]), ks.kernel(ks,0,T,[1 1]), 'half + half = one')
assertRange(ks.kernel(ks,0,T, [1 1]), (exp(-.5) + 1 * 2)/T, 'Test vector');
assertRange(ks.kernel(ks,[0,T], [], [1 1]), (exp(-.5) + 1 * 2)/T, 'Test vec');
assertRange(ks.kernel(ks,[0,1], [], [1 1]), (exp(-.5) * 1.5 + 1.5)/T, 'Test vec');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
ks = kernelFactory('I_exp_int', T);
ksizeList = ks.autoParam(ks, sts)
for cidx = 1:size(ksizeList,1)
    K = computeKernelMatrix(ks, sts, ksizeList{cidx})
end
assertRange(ks.kernel(ks,[],[],1), 1, 'identity exp(0)');
assertRange(ks.kernel(ks,0,0,1), 1, 'identity exp(0)');
assertRange(ks.kernel(ks,1,1,2), 1, 'identity exp(0)');
assertRange(ks.kernel(ks,[1 1.1],[1 1.1],1), 1, 'identity exp(0)');
assertRange(ks.kernel(ks,0,[],1), exp(-T), 'full step');
assertRange(ks.kernel(ks,1,[],1), exp(-(T-1)), 'full step');
assertRange(ks.kernel(ks,[0 1],[],1), exp(-(1+2*2^2)), 'full step');
assertRange(ks.kernel(ks,[0 1],[],5), exp(-(1+2*2^2)/5), 'full step');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
ks = kernelFactory('I_int_exp', T);
ksizeList = ks.autoParam(ks, sts)
for cidx = 1:size(ksizeList,1)
    K = computeKernelMatrix(ks, sts, ksizeList{cidx})
end
assertRange(ks.kernel(ks,[],[],1), 1, 'identity exp(0)');
assertRange(ks.kernel(ks,0,0,1), 1, 'identity exp(0)');
assertRange(ks.kernel(ks,1,1,2), 1, 'identity exp(0)');
assertRange(ks.kernel(ks,[1 1.1],[1 1.1],1), 1, 'identity exp(0)');
assertRange(ks.kernel(ks,0,[],1), exp(-1), 'full step');
assertRange(ks.kernel(ks,1,[],1), 1/3 + exp(-1) * 2 / 3, 'full step');
assertRange(ks.kernel(ks,1,[],2), 1/3 + exp(-1/2) * 2 / 3, 'full step');
assertRange(ks.kernel(ks,[0 1],[],1), exp(-1)/3 + exp(-2^2)*2/3, 'full step');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
ks = kernelFactory('Schrauwen', T, 'Gaussian');
ksizeList = ks.autoParam(ks, sts)
for cidx = 1:size(ksizeList,1)
    K = computeKernelMatrix(ks, sts, ksizeList{cidx})
end
assertRange(ks.kernel(ks,[],[],1), 0, 'empty summation');
assertRange(ks.kernel(ks,[0 1 2 3],[],1), 0, 'empty summation');
assertRange(ks.kernel(ks,0,0,1), 1, 'K(0)');
assertRange(ks.kernel(ks,1,1,1), 1, 'K(0)');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
ks = kernelFactory('spikernel', T, 'Gaussian');
ksizeList = ks.autoParam(ks, sts)
for cidx = 1:size(ksizeList,1)
    K = computeKernelMatrix(ks, sts, ksizeList{cidx})
end
assertRange(norm(K - K'), 1e-30, 'not symmetric!');
EV = eig(K);
assert(min(EV) >= 0, sprintf('not positive definite!! %f', min(EV)));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
ks = kernelFactory('metric-pseudo', T, {'vp', 1}, 'laplacian');

ksizeList = ks.autoParam(ks, sts)
for cidx = 1:size(ksizeList,1)
    K = computeKernelMatrix(ks, sts, ksizeList{cidx})
end

assertRange(ks.kernel(ks,[], [], 1), ks.cmf(0), 'CMF(0) when identical');
assertRange(ks.kernel(ks,[0 0.5], [0 0.5], 1), ks.cmf(0), 'CMF(0) when identical');
assertRange(ks.kernel(ks,[0], [], 1), ks.cmf(1), 'CMF(1) removed');
assertRange(ks.kernel(ks,[0 0.1], [], 1), ks.cmf(2), 'CMF(2) 2 removed');
assertRange(ks.kernel(ks,[0 0.1], [], 2), ks.cmf(1), 'CMF(1) 2 removed scaled');
assertRange(ks.kernel(ks,[0], [0.5], 1), ks.cmf(0.5), 'CMF(1) when difference 1');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
sts = cell(6, 1); % work around of auto kernel size determination problem
sts{1} = [];
sts{2} = [0];
sts{3} = [1];
sts{4} = [1.5];
sts{5} = [1.2];
sts{6} = [0.9];
sts{7} = [1.2 1.6];

ks = kernelFactory('Stratified', T, 'Gaussian');
ksizeList = ks.autoParam(ks, sts)
for cidx = 1:size(ksizeList,1)
    K = computeKernelMatrix(ks, sts, ksizeList{cidx})
end
assertRange(ks.kernel(ks,0,[],1), 0, 'different dimension');
assertRange(ks.kernel(ks,[0 1],[],1), 0, 'different dimension');
assertRange(ks.kernel(ks,[0 1 2],[],1), 0, 'different dimension');
assertRange(ks.kernel(ks,[0 1 2],[1 2],1), 0, 'different dimension');
assertRange(ks.kernel(ks,[],[],1), 1, 'empty identity');
assertRange(ks.kernel(ks,0,0,1), 1, 'K(0)');
assertRange(ks.kernel(ks,1,1,1), 1, 'K(0)');
assertRange(ks.kernel(ks,[0 1],[0 1],3), 1, 'K(0)');
assertRange(ks.kernel(ks,1,2,3), exp(-1/2/3^2), '1D test vector');
assertRange(ks.kernel(ks,[0 1],[1 2],3), exp(-2/2/3^2), '2D test vector');

ks = kernelFactory('Stratified', T, 'Gaussian', 'PoissonScaling');
ksizeList = ks.autoParam(ks, sts)
for cidx = 1:size(ksizeList,1)
    K = computeKernelMatrix(ks, sts, ksizeList{cidx})
end
assertRange(ks.kernel(ks,0,[],1), 0, 'different dimension');
assertRange(ks.kernel(ks,[0 1],[],1), 0, 'different dimension');
assertRange(ks.kernel(ks,[0 1 2],[],1), 0, 'different dimension');
assertRange(ks.kernel(ks,[0 1 2],[1 2],1), 0, 'different dimension');
assertRange(ks.kernel(ks,[],[],1), 1, 'empty identity');
assertRange(ks.kernel(ks,0,0,1), 1, 'K(0)');
assertRange(ks.kernel(ks,1,1,1), 1, 'K(0)');
assertRange(ks.kernel(ks,[0 1],[0 1],3), 1, 'K(0)');
assertRange(ks.kernel(ks,1,2,3), exp(-1/2/(3*ks.scalingFactor(1))^2), '1D test vector');
assertRange(ks.kernel(ks,[0 1],[1 2],3), exp(-2/2/(3*ks.scalingFactor(2))^2), '2D test vector');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
sts = cell(5, 1);
T = 4;
sts{1} = [];
sts{2} = [1];
sts{3} = [2];
sts{4} = [1 2];
sts{5} = [1 2 3];
ks = kernelFactory('Stratified_min', T);
ksizeList = ks.autoParam(ks, sts)
for cidx = 1:size(ksizeList,1)
    K = computeKernelMatrix(ks, sts, ksizeList{cidx})
end
assert(all(eig(K) > 0), 'SPD');
assertRange(ks.kernel(ks,0,[],1), 0, 'different dimension');
assertRange(ks.kernel(ks,[1],[],0), 0, 'different dimension');
assertRange(ks.kernel(ks,[1 2],[],1), 0, 'different dimension');
assertRange(ks.kernel(ks,[1 2],[1 2 2.5],0), 0, 'different dimension');
assertRange(ks.kernel(ks,[],[],1), 1, 'empty identity');
assertRange(ks.kernel(ks,1,1,0), (ks.T-1), 'identical 1D');
assertRange(ks.kernel(ks,[1 2],[1 2],1), (ks.T-1)*(ks.T-2), 'identical 2D');
assertRange(ks.kernel(ks,1,2,0), (ks.T-max(1,2)), '1D test vector');
assertRange(ks.kernel(ks,[1 2],[2 3],0), (ks.T-2) * (ks.T-3), '2D test vector');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
ks = kernelFactory('Stratified_delta', T);
ksizeList = ks.autoParam(ks, sts)
for cidx = 1:size(ksizeList,1)
    K = computeKernelMatrix(ks, sts, ksizeList{cidx})
end
assertRange(ks.kernel(ks,0,[],1), 0, 'different dimension');
assertRange(ks.kernel(ks,[1],[],0), 0, 'different dimension');
assertRange(ks.kernel(ks,[1 2],[],1), 0, 'different dimension');
assertRange(ks.kernel(ks,[1 2],[1 2 2.5],0), 0, 'different dimension');
assertRange(ks.kernel(ks,[],[],1), 1, 'empty identity');
assertRange(ks.kernel(ks,1,1,0), 1, 'same dimension');
assertRange(ks.kernel(ks,[1 2],[1 2],1), 1, 'same dimension');
assertRange(ks.kernel(ks,1,2,0), 1, 'same dimension');
assertRange(ks.kernel(ks,[1 2],[2 3],0), 1, 'same dimension');
assertRange(ks.kernel(ks,[1 2],[3],0), 0, 'different dimension');
