% Generate band limited Gaussian noise

randn('seed', 20110720); rand('seed', 20110720);

nT = 3e5; % Number of time steps
samplingFreq = 1e5;
power = 0.2; % signal power per bin
bias = 3e-2; % bias current per time bin
tr = (1:nT) / samplingFreq; % time range
T = tr(end);

bandLimit = 10;
nSamplesNyquistLimit = bandLimit / samplingFreq * nT;
threshold = 20;

%% First generate white Gaussian noise
x = sqrt(power) * randn(nT, 1) + bias;

subplot(2,1,1); plot(x);
subplot(2,1,2); plot(abs(fftshift(fft(x))));
drawnow;

%% Low-pass filtering with Fourier transform
xf = fft(x);
cutOffBin = ceil(nT / samplingFreq * bandLimit);
xf(cutOffBin:end-cutOffBin+2) = 0;
blx = ifft(xf); % band limited x

subplot(2,1,1); plot(tr, blx);
subplot(2,1,2); plot(abs(fftshift(fft(blx))));
drawnow;

%% Generate spikes
cx = cumsum(blx);
%subplot(2,1,1); plot(tr, cx);
drawnow;
subplot(2,1,2); cla;
currentThreshold = threshold;
st = [];
for t = 1:nT
    if cx(t) >= currentThreshold;
	st(end+1) = t;
	currentThreshold = currentThreshold + threshold;
    end
end
st = st / samplingFreq;
fprintf('Number of spikes generated %d vs Nyquist limit %d\n', length(st), nSamplesNyquistLimit);
st = st(:)';
line([st; st], [0; 1]*ones(size(st)) + [0.1; -0.1]*ones(size(st)), 'Color', 'k');
linkaxes([subplot(2,1,1) subplot(2,1,2)], 'x');
axis tight

subsampleInterval = [1/5 3/5] * T;
nSubsample = 200;
[subsamples, targets, sts, temporalWindow] = timeSeriesRandomSample(blx, st, samplingFreq, nSubsample, subsampleInterval, true);
x = sts;
y = targets;

subplot(2,1,1);
cla
hold all;
plot(tr, blx);
plot(subsamples, targets, 'o');
line(subsamples(1) + [-temporalWindow 0], min(blx) * [1 1], 'LineWidth', 3, 'COlor', [4 4 4]/5);

subsampleInterval = [3/5 5/5] * T;
nSubsample = 200;
[subsamples, targets, sts, temporalWindow] = timeSeriesRandomSample(blx, st, samplingFreq, nSubsample, subsampleInterval, false);
x2 = sts;
y2 = targets;
