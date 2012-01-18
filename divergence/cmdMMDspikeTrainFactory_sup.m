% Estimate the rejection rate of sup-MMD on standard experiments

kernelList = {'wilcoxon', 'count', 'spikernel', 'mci', 'nci1', 'nci2', 'I_exp_int', 'I_int_exp', 'stratified_poisson_scaling', 'reek'}; % chronological ordering
grayColors = (exp((0:11)'/12)-1)/2 * ones(1,3);

% Number of samples to scan for
nList = round(logspace(log10(10), log10(120), 5));
% Number of Monte Carlo for rejection rate
nMC = 100;
alpha = 0.05;

if true % fast run for debugging <==!!
    nMC = 3; nList = [30 50];
    kernelList = {'count', 'mci', 'nci1', 'nci2'};
    verbose = true;
else
    verbose = false;
end 

% rejection or not. actually binary...+NaN
results = nan(length(nList), nMC, length(kernelList));
resultsTC = nan(length(nList), nMC, length(kernelList)); % computation time
mmd2s = nan(length(nList), nMC, length(kernelList));
pvalues = nan(length(nList), nMC, length(kernelList));

% Experiment number
% expNum = 1;
if ~exist('expNum')
    error('First set the experiment number! ex) expNum = 1;');
end

for nIdx = 1:length(nList)
    Nn = nList(nIdx);
    for kMC = 1:nMC
        %% Generate data
        [sts1, sts2] = spikeTrainFactory(expNum, Nn);
        T = sts1.duration;

	%% Do the test
	[pv, res, resTC, mmd2] = divMMDsupBootstrap(sts1(1).data, sts2(1).data, kernelList, T, alpha, verbose);
	pvalues(nIdx, kMC, :) = pv;
	results(nIdx, kMC, :) = res;
	resultsTC(nIdx, kMC, :) = resTC;
	mmd2s(nIdx, kMC, :) = mmd2;
    end % nMC
end % nList

% Save results
save([mfilename '_' datestr(now, 30)]);

%%
resultsBest = squeeze(mean(results, 2)); % output: nIdx x kKernel
nKernels = length(kernelList);

%%
fig = figure(1819); clf; hold on
set(gca, 'FontSize', 10);
bh = bar(resultsBest','edgecolor','none');
set(gca, 'XTick', 1:nKernels);
set(gca, 'FontSize', 6);
set(gca, 'XTickLabel', kernelList);
ylabel('rejection rate');
line([0 nKernels+1], alpha * [1 1], 'LineStyle', '--', 'Color', 'k');
lh = legend([repmat('n = ',length(nList),1), num2str(nList')],'location','NorthWest');
set(lh, 'FontSize', 8);
set(lh, 'DataAspectRatio', [2 1 2])
set(lh, 'box', 'off');

%% save the plot
set(gcf, 'PaperSize', [6 2]);
set(gcf, 'PaperPosition', [0 0 6 2]);
set(gcf, 'PaperUnits', 'Inches');
saveas(gcf, sprintf('power_%d_%s.pdf', expNum, datestr(now,30)));
