% Generate 2nd order Volterra like function from spike trains to real value

randn('seed', 20110723); rand('seed', 20110723);
T = 20;
nSpikes = T * 20; % 20 Hz
st = sort(rand(1, nSpikes)) * T; % Poisson generation
st(diff(st) < 1e-3) = [];
nSpikes = length(st);

%%
% Causal Volterra series
% f(t) = k0 + sum_i k1(t - st(i)) + sum_{i,j} k2(t - st(i), t - st(j))
k0 = 0;
% k1(t) = exp(-t)
k1 = @(t) 200 * -t .* exp(t/0.01);
% k2(t,s) = 0.2 * exp(-(t+s)^2)
k2 = @(t,s) (-0.2 * exp(-(t+s).^2/0.01));
% Compute value at t given spike timings in the past
f = @(st1, stst1, stst2) (k0 + sum(k1(st1)) + sum(sum(k2(stst1, stst2))));
tr = (1:1000)/1000*T;

clf;
subplot(1,2,1);
[stst1, stst2] = meshgrid(linspace(0, -.5));
mesh(-stst1, -stst2, k2(stst1, stst2));
xlabel('t_1'); ylabel('t_2'); zlabel('k_2(t_1, t_2)');
subplot(1,2,2);
plot(linspace(0, .5, 1000), k1(linspace(0, -.5, 1000)));
xlim([0 .2]); set(gca, 'Ytick', [])
xlabel('t'); ylabel('k_1(t)')

%%
[fval fval1 fval2] = volterraEval(k0, k1, k2, st, tr);

%%
clf;
subplot(2,1,1); hold all;
plot(tr, fval); plot(tr, fval1); plot(tr, fval2);

subplot(2,1,2);
cla; plotRaster({st}, 'k');

%%
samplingFreq = 1/(tr(2) - tr(1));
subsampleInterval = [1/5 4/5] * T;
nSubsample = 200;
temporalWindow = 0.2;
[subsamples1, targets, sts, temporalWindow] = timeSeriesRandomSample(fval, st, samplingFreq, nSubsample, subsampleInterval, true, temporalWindow);
x = sts;
targets = volterraEval(k0, k1, k2, st, subsamples1);
y = targets;

subplot(2,1,1);
cla
hold all;
plot(tr, fval);
plot(subsamples1, targets, 'o');
line(subsamples1(1) + [-temporalWindow 0], min(fval) * [1 1], 'LineWidth', 3, 'COlor', [4 4 4]/5);

subsampleInterval = [4/5 5/5] * T;
nSubsample = 200;
[subsamples, targets, sts, temporalWindow] = timeSeriesRandomSample(fval, st, samplingFreq, nSubsample, subsampleInterval, false, temporalWindow);
targets = volterraEval(k0, k1, k2, st, subsamples);
x2 = sts;
y2 = targets;

linkaxes([subplot(2,1,1) subplot(2,1,2)], 'x');