function [x binCenters binEdges] = binSpikeTrains(spikeTrains, T, binSize)
% x = binSpikeTrains(spikeTrains, T, binSize)
% Bin a set of spike trains represented as a set of timings.
% It assumes the spike timings are between 0 and T.
%
% Output
%   x: (NxM) M spike trains, N bins
%   binCenters: (Nx1) bin centers
%   binEdges: ((N+1)x1) bin edges
%
% Copyright 2011-2012 Memming. All rights reserved

N = ceil(T / binSize);
M = length(spikeTrains);

x = zeros(N, M);
binEdges = 0:binSize:T;
binCenters = (binEdges(1:end-1) + binEdges(2:end))/2;

for kk = 1:M
    h = histc(spikeTrains{kk}, binEdges);
    if isempty(h)
        continue;
    end

    if N == length(binEdges)
	x(:, kk) = h;
    else
	x(:, kk) = h(1:end-1);
	x(end, kk) = x(end, kk) + h(end); % include the spikes on the boundary in the last bin
    end
end
