function [sts, targets, T, targetRange] = regressionFactory(expNum, N)
% Create an artificial regression problem
% Input
%   expNum: experiment number
%   N: number of samples to draw
% Output
%   sts: {Nx1} spike trains
%   targets: (Nx1) target values
%   T: length of spike trains
%   targetRange: (Dx1) range of target values
%
% $Id$

% We need a set of regression problems from spike train to real value.
% These are decoding tasks in neuroscience, so simply we can design an
% ecoding model and reverse it.
% 
% Find $f$ such that $y = f(x)$.
% 
% Experiment 1: simple rate coding
% 
% x \in [0 40]
% y ~ PoissionProcess(x * [0 T])
% 
% Experiment 2: rate coding with tuning curve
% 
% x \in [0 2\pi]
% bias = 20
% gain = 15
% y ~ PoissionProcess(gain * cos(x)+bias * [0 T])
% 
% Experiment 3: ISI coding with Gamma shape parameter
% 
% x \in [0 10]
% y ~ RenewalProcess(Gamma(x, lambda))
% 
% Experiment 4: temporally modulated rate
% 
% x \in [0 2 \pi]
% t \in [0 2 \pi]
% y ~ PoissonProcess(cos((t + x)/T))
%
% Experiment 8: same count distribution, random locations
% - Original signal/noise demo
%
% Experiment 9: perfectly periodic spike trains with random phase
% -> This tests perfect noise/signal seperation 

sts = cell(N, 1);
T = 1;

switch expNum
case 1
    targetRange = [0 40];
    targets = targetRange(1) + rand(N, 1) * diff(targetRange);
    for k = 1:N
	sts{k} = sort(rand(poissrnd(targets(k)),1)) * T;
    end
case 1.5
    targetRange = [0 40];
    targets = ceil(targetRange(1) + rand(N, 1) * diff(targetRange));
    for k = 1:N
	sts{k} = sort(rand(targets(k), 1)) * T;
    end
case 2
    targetRange = [0 2*pi];
    targets = rand(N, 1) * targetRange(2);
    bias = 20;
    gain = 15;
    assert(bias >= gain, 'bias must be larger than gain');
    for k = 1:N
	st = sort(rand(poissrnd(bias * T),1));
	sts{k} = cosineTimeRescaling(st, gain, bias, phase, T);
    end
case 5
    % Generates spike at a given location with small jitter
    % Designed to be a simple regression problem
    jitterSD = 0.1;
    targets = rand(N,1);
    targetRange = [0 1];
    for k = 1:N
	sts{k} = targets(k) + jitterSD*randn;
	if sts{k} <= 0
	    sts{k} = [];
	end
    end
case 6
    % Extract the location of the first spike with certain jitter
    % Designed to perform regression under redundancy
    minSpike = 3;
    maxSpike = 7;
    jitterSD = 0.01;
    targets = zeros(N,1);
    targetRange = [-Inf Inf]; % <-- normal distribution
    for k = 1:N
	sts{k} = sort(rand(minSpike + ceil((maxSpike+1-minSpike)*rand), 1));
	targets(k) = min(sts{k}) + jitterSD*randn;
    end

case 7
    % Serially correlated process
    % Predict the correlation
    jitterSD = 0.01;
    targets = rand(N,1);
    targetRange = [0 1];
    for k = 1:N
	sts{k}(1) = 0.5 +  jitterSD * randn;
	sts{k}(2) = targets(k) * sts{k}(1) + 1;
    end

case 8
    targetRange = linspace(0, T, 20);
    dTarget = targetRange(2) - targetRange(1);
    countDistribution = [0 15 30 5 10 20 7 5 3 1 1];
    countDistribution = countDistribution(:) / sum(countDistribution);
    countRange = 0:(length(countDistribution)-1);
    NN = randsample(countRange, N,true,countDistribution);

    targets = zeros(N,1);
    for k = 1:N
	targets(k) = targetRange(ceil(rand * (length(targetRange)-1)));
	sts{k} = targets(k) + sort(rand(NN(k),1)) * dTarget;
    end
    targetRange = [targetRange(1) targetRange(end)];

case 9
    targets = T/10 + T * rand(N,1)/20; % period
    targetRange = T/10 + [0 1/20];
    for k = 1:N
       	period = targets(k);
	phi = rand * period; % The phase
	sts{k} = phi:period:T;
    end
otherwise
    error('No such experiment number [%d]', expNum);
end

end % function

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [st] = cosineTimeRescaling(st, gain, bias, phase, T)
% Use the time rescaling theorem to transform a Poisson

% t \in [0 2 \pi]
% \lambda(t) = gain * cos(t + phase) + bias
% \Lambda(t) = gain * sin(t + phase) + bias * t

end
