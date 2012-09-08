function [metric] = metricFactory(metricParam)
% [metric] = metricFactory(metricParam);
% Returns a spike train metric function handle @(st1, st2) -> non-neg real
%
% Input:
%    metricParam: cell array of paramters. 
%	First is the name of the distance: VP or DSM
%	The rest are usually real valued parameters. q for VP, p and q for DSM
%
% Copyright 2011 Memming. All rights reserved.

if nargin < 1 || isempty(metricParam)
    metricParam = {'vp', 1}; % Default
end

if isstr(metricParam)
    metricName = metricParam;
elseif iscell(metricParam)
    metricName = metricParam{1};
else
    error('Metric factory takes a cell with first element being the name');
end

switch lower(metricName)
case {'vp', 'victor-purpura'}
    q = metricParam{2};
    metric = @(st1, st2) dsm_metric(st1, st2, 1, q);
case {'dsm', 'lp', 'lpalignment'}
    q = metricParam{2};
    p = metricParam{3};
    metric = @(st1, st2) dsm_metric(st1, st2, p, q);
otherwise
    warning('Unknown metric name');
end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function d = dsm_metric(X,Y,p,q)
% Input:
%   X: spike timings
%   Y: spike timings
%   p: L-p distance (p = 1 is the Victor-Purpura metric)
%   q: cost of moving a spike by a unit distance in time
%
% Ref: Alexander J. Dubbs, Brad A. Seiler and Marcelo O. Magnasco,
% A Fast Lp Spike Alignment Metric, Neural Computation 22 (2010): 2785--2808
% Implementation by Alexander Dubbs (personal communication) 
% copyright status unknown

M = initM(X,Y,p,q);
v = 0;
for j = 1:size(M,2)
    shifted = 1;
    while shifted == 1
        [M{j},shifted,length_shift] = shiftM(M{j},X,Y,p,q);
        v = v - length_shift;
    end
end
d = (v + length(X) + length(Y))^(1/p);

end

function M = initM(X,Y,p,q)
M = cell(0);
xcount = 0;
ycount = 0;
Mcount = 1;
M{1} = [];
thresh = 2^(1/p)/q;
lastval = NaN;
while ~and(isempty(X),isempty(Y))
    if isempty(X)
        if Y(1)-lastval > thresh
            Mcount = Mcount + 1;
            M{Mcount} = [];
        end
        ycount = ycount + 1;
        M{Mcount} = [M{Mcount},[NaN;ycount]];
        lastval = Y(1);
        Y(1) = [];
    else
        if isempty(Y)
            if X(1)-lastval > thresh
                    Mcount = Mcount + 1;
                    M{Mcount} = [];
            end
            xcount = xcount + 1;
            M{Mcount} = [M{Mcount},[xcount;NaN]];
            lastval = X(1);
            X(1) = [];
        else
            if X(1) <= Y(1)
                if X(1)-lastval > thresh
                    Mcount = Mcount + 1;
                    M{Mcount} = [];
                end
                    xcount = xcount + 1;
                    M{Mcount} = [M{Mcount},[xcount;NaN]];
                    lastval = X(1);
                    X(1) = [];
            else
                if Y(1)-lastval > thresh
                    Mcount = Mcount + 1;
                    M{Mcount} = [];
                end
                ycount = ycount + 1;
                M{Mcount} = [M{Mcount},[NaN;ycount]];
                lastval = Y(1);
                Y(1) = [];
            end
        end
    end
end

end

function [M_new,shifted,length_max] = shiftM(M,X,Y,p,q)
q_p = q^p;
length_max = 0; % NB If a shift is of nonpositive length, it will never be made.
length_counter = 0;
up_down_old = 2; % 1 denotes last NaN was up, 0 denotes it was down,
                 % 2 means we haven't encountered a NaN yet.
shifted = 0;
Left_old = [];
Left = [];
Right = [];
for i = 1:size(M,2)
    if ~or(isnan(M(1,i)),isnan(M(2,i)))
        if up_down_old == 0
            length_counter = length_counter - q_p*abs(X(M(1,i-1))-Y(M(2,i)))^p + q_p*abs(X(M(1,i))-Y(M(2,i)))^p;
        else
            if up_down_old == 1
                length_counter = length_counter - q_p*abs(X(M(1,i))-Y(M(2,i-1)))^p + q_p*abs(X(M(1,i))-Y(M(2,i)))^p;
            end
        end
    else
        if isnan(M(1,i))
            if up_down_old == 0
                length_counter = length_counter - q_p*abs(X(M(1,i-1))-Y(M(2,i)))^p + 2;
                if length_counter > length_max
                    length_max = length_counter;
                    Left = Left_old;
                    Right = i;
                end
            end
            length_counter = 0;
            Left_old = i;
            up_down_old = 1;
        else
            if up_down_old == 1
                length_counter = length_counter - q_p*(X(M(1,i))-Y(M(2,i-1)))^p + 2;
                if length_counter > length_max
                    length_max = length_counter;
                    Left = Left_old;
                    Right = i;
                end
            end
            length_counter = 0;
            Left_old = i;
            up_down_old = 0;
        end
    end
end
if isempty(Left)
    M_new = M;
else
    shifted = 1;
    row1 = M(1,:);
    row2 = M(2,:);
    if isnan(row1(1,Left))
        row1(Left) = [];
        row2(Right) = [];
    else
        row1(Right) = [];
        row2(Left) = [];
    end
    M_new = [row1;row2];
end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
