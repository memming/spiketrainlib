function sts = genNegBinHomogeneous(meanCount, varCount, T, N)
% Count distribution follows negative binomial.
% Negative binomial with two parameters r and p
% mean = p*r/(1-p)
% var = p*r/(1-p)^2
% Note that MATLAB stat toolbox's nbinrnd takes r and (1-p)

p = 1 - meanCount / varCount;
if p <= 0
    error('Variance is too low for negative binomial (FF > 1)');
end
r = meanCount / p * (1 - p);
cdrnd = @() nbinrnd(r, 1 - p);

sts = cell(N, 1);
for k = 1:N
    sts{k} = sort(rand(cdrnd(), 1)) * T;
end
