function val = virspike(ks, s1, s2, param)

s1 = s1 / ks.T; s2 = s2 / ks.T;

if ~isempty(s1)
    t1 = zeros(2*length(s1)+1,1);
    t1(end) = 1; t1(2:2:end-1) = s1;
    t1(3:2:end-1) = 0.5*(s1(1:end-1) + s1(2:end));
    a1 = zeros(length(t1),1);
    a1(2:2:end-1) = 0;
    a1(3:2:end-1) = 0.5*(s1(2:end) - s1(1:end-1));
    a1(1) = s1(1); a1(end) = 1 - s1(end);
else
    t1(1,1) = 0; t1(2,1) = 1;
    a1(1,1) = 1; a1(2,1) = 1;
end
if ~isempty(s2)
    t2 = zeros(2*length(s2)+1,1);
    t2(end) = 1; t2(2:2:end-1) = s2;
    t2(3:2:end-1) = 0.5*(s2(1:end-1) + s2(2:end));
    a2 = zeros(length(t2),1);
    a2(2:2:end-1) = 0;
    a2(3:2:end-1) = 0.5*(s2(2:end) - s2(1:end-1));
    a2(1) = s2(1); a2(end) = 1 - s2(end);
else
    t2(1,1) = 0; t2(2,1) = 1;
    a2(1,1) = 1; a2(2,1) = 1;
end

t = sort(unique([t1;t2]));
a = interp1(t1,a1,t);
b = interp1(t2,a2,t);

m1 = (a(2:end) - a(1:end-1)) ./ (t(2:end) - t(1:end-1));
m2 = (b(2:end) - b(1:end-1)) ./ (t(2:end) - t(1:end-1));
b1 = a(1:end-1); b2 = b(1:end-1);

val = sum((m1.*m2).*(t(2:end) - t(1:end-1)).^3/3 ...
        + (m1.*b2 + m2.*b1).*(t(2:end) - t(1:end-1)).^2/2 ...
        + (b1.*b2).*(t(2:end) - t(1:end-1)));
