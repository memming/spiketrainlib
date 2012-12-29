% KPCA for divergecne (MMD) is not the best thing to see
% what "orthogonal features" MMD is most sensitive to.
% It would be better to compute MMD for each PC direction
% and display the ones with highests difference.

expNum = 3;
kernelList = {'mci', 'nci1', 'spikernel'}; %'stratified Poisson scaling';

N = 6 * 32; % number of samples per class (total data = 2*N)
testSize = 0.05;

[sts1, sts2, desc] = spikeTrainFactory(expNum, N);
T = sts1.duration;
idx1 = 1:N';
idx2 = (N+1):(2*N)';
allsts = cell(2*N, 1);
allsts(idx1) = sts1(1).data(:);
allsts(idx2) = sts2(1).data(:);
sts = allsts([idx1(:); idx2(:)]);

%% plot spike trains
figRaster = figure;
subplot(2,2,1); plotRaster(sts1.data, [0 0 0], [], false); title(desc);
subplot(2,2,3); plotRaster(sts2.data, [0 0 0], [], false); title(desc);
subplot(2,2,2); plotRaster(sts1.data, [0 0 0], [], true); title('sorted');
subplot(2,2,4); plotRaster(sts2.data, [0 0 0], [], true); title('sorted');

%% 
figPCA = figure;
nKernels = length(kernelList);

for kKernel = 1:nKernels
kernelName = kernelList{kKernel};
ks = kernelFactory(kernelName, T);
fprintf('Kernel [%s]\n', ks.name);

if false
    %% just 1 kernel size
    ksizeList = ks.autoParam(ks, sts);
    cidx = ceil(length(ksizeList)/2);
    KM = computeKernelMatrix(ks, sts, ksizeList{cidx});
    [rejected, mmd2, threshold, d2] = mmdBootstrap(KM, idx1, idx2, testSize);
else
    %% do it for all kernel sizes and choose the best one
    ksizeList = ks.autoParam(ks, sts);
    nks = length(ksizeList);
    fprintf('scanning %d kernel sizes\n', nks);

    rs = rand('seed'); % <-- !!!
    KM = zeros(2*N, 2*N, nks);
    mmd2All = zeros(nks, 1);
    d2All = zeros(nks, 9999); % TODO magic number in mmd_bootstrapping
    for ksidx = 1:nks
	%% KPCA and sort eigenvalues
	KM(:, :, ksidx) = computeKernelMatrix(ks, sts, ksizeList{ksidx});
	fprintf('%d/%d\r', ksidx, nks);

	rand('seed', rs); % hack to get the same bootstrap samples (TODO)
	[rejected, mmd2_instance, threshold, d2_instance] = ...
	    mmdBootstrap(KM(:, :, ksidx), idx1, idx2, testSize);
	mmd2All(ksidx) = mmd2_instance;
	d2All(ksidx, :) = d2_instance;
    end
    %% Compute the sup MMD value
    [mmd2 midx] = max(mmd2All); d2 = max(d2All, [], 1);
    fprintf('\nkernel size chosen [%d]\n', midx);
    disp(ksizeList{midx});
    KM = KM(:, :, midx);
end % kernel size if

p1 = sum(d2 > mmd2) / length(d2);
p2 = 1 - sum(d2 < mmd2) / length(d2);
p = (p1 + p2) / 2;

%% max dimension (dividing by small eigVal is not numerically stable)
md = min(N*2, 100);
[proj, eigVal, eigVec] = KPCA(KM, md); % sorts internally
proj = bsxfun(@times, proj(:, 1:md), 1 ./ sqrt(eigVal(1:md))'); % truncate

%% Compute component wise MMD
d2 = ((mean(proj(idx1, 1:md)) - mean(proj(idx2, 1:md)))'.^2);
mmd2alt = sum(d2);

%% Clean up numerical instability
eigVal(abs(eigVal) < 5 * eps) = 0;
eigVal = real(eigVal); % remove imaginary values
d2 = real(d2);
mmd2alt = real(mmd2alt);
fprintf('MMD2 %g [p=%f] MMD2 alt (via KPCA) %g\n', mmd2, p, mmd2alt);

%% plot MMD2 components
subplot(nKernels,2,1 + (kKernel-1)*2); cla; hold all
plot(eigVal / sum(eigVal), 'or', 'MarkerSize', 8)
plot(d2 / mmd2alt, '.k', 'MarkerSize', 8)
xlabel('sorted PC');
ylabel('normalized MMD2');
xlim([0 20]);
title('Normalized eigenvalues of P+Q and MMD2 per dimension');
legend('spectrum', 'MMD2 component');

%% Plot KPCA for the highest PCs
[~, sidx] = sort(d2, 'descend');
subplot(nKernels,2,2 + (kKernel-1)*2); cla; hold on;
plot(proj(idx1, sidx(1)), proj(idx1, sidx(2)), 'dk');
plot(mean(proj(idx1, sidx(1))), mean(proj(idx1, sidx(2))), 'dk', 'MarkerFaceColor', 'k');
plot(proj(idx2, sidx(1)), proj(idx2, sidx(2)), 'or');
plot(mean(proj(idx2, sidx(1))), mean(proj(idx2, sidx(2))), 'sr', 'MarkerFaceColor', 'r');
xlabel(sprintf('PC [%d]', sidx(1)));
ylabel(sprintf('PC [%d]', sidx(2)));
title(desc);
axis tight
drawnow

end % kKernel
