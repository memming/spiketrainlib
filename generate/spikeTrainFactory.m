function [spikeTrains1, spikeTrains2, desc] = spikeTrainFactory(expNum, N)
% TODO rename this

M = 1;
switch(expNum)
case 1
desc = 'Two correlated vs independent APs';
spikeTrains1 = genTwoAPex(N, M, struct('jitter',10e-3,'type','correlated'));
spikeTrains2 = genTwoAPex(N, M, struct('jitter',10e-3,'type','uncorrelated'));
case 2
desc = 'Homogeneous Poisson processes';
fprintf('=== Experiment 2: %s\n', desc);
homogeneousParam1.tOffset = 0.25;
homogeneousParam1.duration = 0.50;
homogeneousParam1.lambda = 3;
homogeneousParam2 = homogeneousParam1;
homogeneousParam2.lambda = 5;
spikeTrains1 = genHomogeneousPoisson(N, M, homogeneousParam1);
spikeTrains2 = genHomogeneousPoisson(N, M, homogeneousParam2);
case 3
desc = 'PTST vs equi-Poisson process';
fprintf('=== Experiment 3: %s\n', desc);
PTSTParams.L = 3;
PTSTParams.T = 1;
PTSTParams.type = 'PTST';
[spikeTrains1, PoissonPTSTParams] = genPTST(N, M, PTSTParams);
PoissonPTSTParams.type = 'equPoisson'
spikeTrains2 = genPTST(N, M, PoissonPTSTParams);
case 4
desc = 'Renewal vs Serially correlated';
fprintf('=== Experiment 4: %s\n', desc);
renewalParam = struct('T',125e-3,'mISI',50e-3,'urISI',5e-3,'type','renewal');
serialCorrParam = struct('T',125e-3,'mISI',50e-3,'urISI',5e-3,'type','correlated');
spikeTrains1 = genSerialCorr(N, M, renewalParam);
spikeTrains2 = genSerialCorr(N, M, serialCorrParam);
case 5
desc = 'Stratified specification (Gaussian)';
fprintf('=== Experiment 5: %s\n', desc);
stratParam1.D = 2;
stratParam1.p = [0.5 0.5];
stratParam1.T = 1;
stratParam1.mu{1} = 0.6;
stratParam1.sigma{1} = 20e-3^2;
stratParam1.mu{2} = [0.4 0.6];
stratParam1.sigma{2} = [20e-3^2, 0; 0, 20e-3^2];
stratParam2 = stratParam1;
stratParam2.mu{2} = [0.4 0.5];
spikeTrains1 = genStratifiedGaussian(N, M, stratParam1);
spikeTrains2 = genStratifiedGaussian(N, M, stratParam2);
case 6
desc = 'Step Poisson process with different rate';
fprintf('=== Experiment 6: %s\n', desc);
stepPoissonParam1 = struct('sectionLength', 100e-3, 'lambda1', 3, 'lambda2', 2);
stepPoissonParam2 = struct('sectionLength', 100e-3, 'lambda1', 2, 'lambda2', 3);
spikeTrains1 = genStepPoisson(N, M, stepPoissonParam1);
spikeTrains2 = genStepPoisson(N, M, stepPoissonParam2);
case 7
desc = 'Renewal process with different interval structure';
fprintf('=== Experiment 7: %s\n', desc);
renewalParam1.T = 1;
renewalParam1.rate = 10;
renewalParam2 = renewalParam1;
renewalParam1.shape = 3;
renewalParam2.shape = 0.5;
spikeTrains1 = genRenewalGamma(N, M, renewalParam1);
spikeTrains2 = genRenewalGamma(N, M, renewalParam2);
case 8
desc = 'Renewal process vs Poisson';
fprintf('=== Experiment 8: %s\n', desc);
renewalParam1.T = 1;
renewalParam1.rate = 10;
%renewalParam1.shape = 2;
renewalParam1.shape = 10;
spikeTrains1 = genRenewalGamma(N, M, renewalParam1);
homogeneousParam1.tOffset = 0;
homogeneousParam1.duration = 1;
homogeneousParam1.lambda = 10;
spikeTrains2 = genHomogeneousPoisson(N, M, homogeneousParam1);
case 9
desc = 'PTPP with high noise ';
fprintf('=== Experiment 9: %s\n', desc);
PTSTParams.L = 3;
PTSTParams.T = 3;
PTSTParams.mu = [0.5 1.5 2.5]';
PTSTParams.sigma = 0.08 * ones(3, 1);
PTSTParams.p = 0.8 * ones(3, 1);
PTSTParams.type = 'PTST';
PTSTParams.lambda = 4; % <-- POISSON NOISE LEVEL
spikeTrains1 = genPTST(N, M, PTSTParams);
PTSTParams.L = 2;
PTSTParams.T = 3;
PTSTParams.sigma = 0.08 * ones(2, 1);
PTSTParams.p = 0.8 * ones(2, 1);
PTSTParams.mu = [1 2]'; % <-- DIFFERENT POSITIONS
PTSTParams.lambda = 9; % <-- POISSON NOISE LEVEL
PTSTParams.type = 'equPoisson';
spikeTrains2 = genPTST(N, M, PTSTParams);
case 10
desc = 'periodic spike trains (PHASE difference)';
fprintf('=== Experiment %d: %s\n', expNum, desc);
N2 = N * 2;
param.T = 1;
param.countList = 3:8; %linspace(0.15, 0.3, param.nPeriod);
param.nPhase = N2 / length(param.countList);
[sts phaseIdx periodIdx] = genPeriodic(param, N2);
spikeTrains1(1).duration = param.T;
spikeTrains1(1).data = sts(phaseIdx <= param.nPhase/2);
spikeTrains2(1).duration = param.T;
spikeTrains2(1).data = sts(phaseIdx > param.nPhase/2);
case 11
desc = 'periodic spike trains (PERIOD difference)';
fprintf('=== Experiment %d: %s\n', expNum, desc);
N2 = N * 2;
param.T = 1;
%param.periodList = linspace(0.15, 0.3, param.nPeriod);
%param.phaseList = linspace(0, 1, param.nPhase);
param.countList = 3:8; %linspace(0.15, 0.3, param.nPeriod);
param.nPhase = N2 / length(param.countList);
[sts phaseIdx periodIdx] = genPeriodic(param, N2);
spikeTrains1(1).duration = param.T;
spikeTrains1(1).data = sts(periodIdx <= length(param.countList)/2);
spikeTrains2(1).duration = param.T;
spikeTrains2(1).data = sts(periodIdx > length(param.countList)/2);
%{
param.nPeriod = round(sqrt(N));
param.nPhase = round(N / param.nPeriod);
param.T = 1;
param.periodList = linspace(0.1, 0.15, param.nPeriod);
param.phaseList = linspace(0, 1, param.nPhase);
spikeTrains1(1).duration = param.T;
spikeTrains1(1).data = genPeriodic(param, N);
param.periodList = linspace(0.15, 0.3, param.nPeriod);
spikeTrains2(1).duration = param.T;
spikeTrains2(1).data = genPeriodic(param, N);
%}
otherwise
    error('No such experiment number');
end
fprintf('Data generation complete...[N = %d, M = %d]\n', N, M);
