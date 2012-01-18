function [fval fval1 fval2] = volterraEval(k0, k1, k2, st, tr)

fval1 = zeros(length(tr),1);
fval2 = zeros(length(tr),1);
for tidx = 1:length(tr);
    t = tr(tidx);
    st1 = st(st < t) - t;
    [stst1, stst2] = meshgrid(st1);
    fval1(tidx) = sum(k1(st1));
    k2val = k2(stst1, stst2);
    fval2(tidx) = sum(k2val(:));
end
fval = k0 + fval1 + fval2;
fprintf('%f, %f\n', var(fval1), var(fval2));