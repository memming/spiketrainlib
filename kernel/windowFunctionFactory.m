function [wfh dwfh] = windowFunctionFactory(T, windowFunctionName)
% Returns a window function f
% Typically, f(0) = f(T) = 0, and f(x) >= 0 for x in [0 T]
%
% $Id: windowFunctionFactory.m 92 2011-09-07 22:36:23Z memming $
% Copyright 2011 Memming. All rights reserved.

if isempty(windowFunctionName)
    warning('Window function is not given.');
    windowFunctionName = 'quadratic';
end

switch lower(windowFunctionName)
    case {'quadratic'}
	wfh = @(t) (-t .* (T - t));
	dwfh = @(t) (2*t - T);
    case {'sin'}
	wfh = @(t) (sin(pi*t/T));
	dwfh = @(t) (pi/T*cos(pi*t/T));
    case {'constant', '1', 'rectangular'}
	wfh = @(t) ones(size(t));
	dwfh = @(t) zeros(size(t));
    case {'inverse gaussian', 'gaussian bump'}
	% http://en.wikipedia.org/wiki/Bump_function
	wfh = @(t) exp(-1./(1-(2*t/T - 1).^2));
	dwfh = @(t) T^2 * exp(-1./(1-(2*t/T - 1).^2)) .* (T - 2*t) ./ (4*t.^2*(T-t).^2);
    otherwise
        error('Unknown window function name [%s]', windowFunctionName);
end % switch
