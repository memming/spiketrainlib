%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fprintf('||||  expNum = 3\n');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
expNum = 3;
param.T = 10;
param.N = 100;
param.lambda = [3 4];
[sts, feature] = singleFeaturePointProcessFactory(expNum, param);

N1 = cellfun('length', sts(:, 1)); N2 = cellfun('length', sts(:, 2));
Nrange = 0:max(max(N1), max(N2));
fprintf('Mean count: (%f, %f) %f, %f\n', param.lambda(1), param.lambda(2), mean(N1), mean(N2));
fprintf('Var count: (%f, %f) %f, %f\n', param.lambda(1), param.lambda(2), var(N1), var(N2));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%      expNum = 1, 2                                %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
param.T = 2;
param.N = 50;
param.meanCount = 4;
param.varCount = [5 10];

for expNum = 1:2
fprintf('||||  expNum = %d\n', expNum);
[sts, feature] = singleFeaturePointProcessFactory(expNum, param);

N1 = cellfun('length', sts(:, 1)); N2 = cellfun('length', sts(:, 2));
Nrange = 0:max(max(N1), max(N2));
fprintf('Mean count: (%f) %f %f\n', param.meanCount, mean(N1), mean(N2));
fprintf('Var count: (%f %f) %f %f\n', param.varCount(1), param.varCount(2), var(N1), var(N2));

figure(131); clf;
subplot(2,2,1);
plotRaster(sts(:, 1));
subplot(2,2,2);
hist(N1, Nrange);
xlim([-1 Nrange(end)+1]);

subplot(2,2,3);
plotRaster(sts(:, 2));
subplot(2,2,4);
hist(N2, Nrange);
xlim([-1 Nrange(end)+1]);
drawnow; pause(0.5); beep;
end

% divMMDsupBootstrap(sts(:, 1), sts(:, 2), {'count', 'spikernel', 'stratified poisson scaling', 'nCI1'}, param.T, 0.1, true);


beep
