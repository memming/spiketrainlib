function A = meanStkernelCount(hyp, x, i)
% Mean function that is linearly dependent on the total count.
% Interface for GPML library (v3.1)

if nargin < 2, A = '1'; return; end

M = cellfun('length', x);
M = M(:);

if nargin < 3
    A = M * hyp;
else
    A = M; % derivative
end
