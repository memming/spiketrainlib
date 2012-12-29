% use bandlimitedGaussian

subsampleInterval = [0 T];
nSubsample = 500;
temporalWindow = []; %0.05;
[subsamples, targets, sts, temporalWindow] = timeSeriesRandomSample(blx, st, samplingFreq, nSubsample, subsampleInterval, false, temporalWindow);

switch 1
case 1
    targets = zscore(targets);
    options.ks = kernelFactory('schoenberg', temporalWindow, 'gaussian');
    options.ksize = [0.010 0.5];
    options.learningRate = 0.5;
case 2
    options.ks = kernelFactory('mCI', temporalWindow);
    options.ksize = 0.005;
    options.learningRate = 0.001;
end
%ksizes = options.ks.autoParam(options.ks, sts);
%options.ksize = ksizes{round(numel(ksizes)/2)};

state = stklms(options);
tic;
KM = computeKernelMatrix(state.ks, sts, state.ksize); % precompute kernel matrix
toc

%% figure setup
figure(6161); clf; 
subplot(5,1,5);
line([st; st], [0; 1]*ones(size(st)), 'Color', 'k');
set(gca, 'XTick', [], 'YTick', [], 'box', 'off');
subplot(5,1,1:4); hold on;
set(gca, 'TickDir', 'out');

%% Fast training with precomputed matrix
state.coeff = state.learningRate * (targets(1) - KM(1,1));
state.x = sts;
ypred = zeros(nSubsample, 1);
ypred1 = ypred;
for k = 2:nSubsample
    yhat = sum(state.coeff .* KM(k, 1:(k-1)));
    ypred1(k) = yhat;
    err = targets(k) - yhat;
    state.coeff(end+1) = state.learningRate * err;
    state.n = state.n + 1;

     for kk = 1:nSubsample
 	ypred(kk) = sum(state.coeff .* KM(kk, 1:k));
     end
     cla; hold on;
     plot(subsamples, targets, 'k');
     plot(subsamples, ypred, 'r');
     plot(subsamples(k), targets(k), 'ok');
     drawnow
end
plot(subsamples, ypred1, 'b');

set(gcf, 'PaperUnits', 'inches', 'PaperSize', [12 4], 'PaperPosition', [0 0 12 4]);
saveas(gcf, 'KLMS_example.pdf');

%% This is SLOW (kernel computed again every time)
% ypred = zeros(nSubsample, 1);
% for k = 1:nSubsample
%     [state, yhat] = stklms(state, sts{k}, targets(k));
% 
%     if k > 100
%     for kk = 1:nSubsample
% 	ypred(kk) = stklms(state, sts{kk});
%     end
%     cla; hold on;
%     plot(subsamples, targets, 'k');
%     plot(subsamples, ypred, 'r');
%     plot(subsamples(k), targets(k), 'ok');
%     beep; pause;
%     end
% end
