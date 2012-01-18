function candidates = autoKernelSizeScalar(scaledQuantityHandle, spikeTrains)
% This only works in the case of scalar
% Input:
%   scaledQuantityHandle: function handle for computing the quantity that needs
%			to be scaled by a linear kernel size
%   spikeTrains: {Nx1} set of spike trains that the kernel will be applied for
% Output:
%   candidates: (3x1) suggested kernel sizes (strictly positive)
%
% Copyright 2010-2012 Memming. All rights reserved.

%% Subsampling at most 100 spike train pairs
nSpikeTrains = numel(spikeTrains);
M = min(200, nSpikeTrains);
k1 = randperm(nSpikeTrains); k1 = k1(1:M);
k2 = randperm(nSpikeTrains); k2 = k2(1:M);

%% Computing the quantity for the subsampled pairs
d = [];
for k = 1:M
    dd = scaledQuantityHandle(spikeTrains{k1(k)}, spikeTrains{k2(k)});
    d = [d; dd(:)];
end

%% Choose the quantiles
if all(d == 0)
    warning('All pairwise scaled quantity is zero!');
    candidates = 0; % TODO Change 0 to NaN
    return;
else
    d = d(d~=0);
    candidates = quantile(d, [0.1; 0.5; 0.9]); % Strictly positive
    candidates = [ candidates(1)/2; candidates; candidates(end) * 2 ];
    %candidates = quantile(d, [0.1 0.5 0.9]); % Strictly positive
    %candidates = candidates(:);
    candidates = median(d) * 2.^[-2:2]'; % Median
    candidates = unique(candidates(:));
    % TODO replace above with below
    % OR better yet, why don't we just leave the duplicates and leave it
    % for the computation part to take care of it?
    %{
    bidx = true(size(candidates));
    [~, uidx] = unique(candidates(:), 'first');
    bidx(uidx) = false;
    candidates(bidx) = NaN;
    %}
end

candidates = mat2cell(candidates, ones(size(candidates,1),1), size(candidates,2));
