function [out, yhat] = stklms(state, x, y)
% Spike train to real value function.
% Initialization mode: state = stklms(options);
% Training mode:       [state, yhat] = stklms(state, x, y);
% Prediction mode:     yhat = stklms(state, x);

if nargin < 1, error('At least one argument is required'), end
if nargin < 2 % initialize
    if ~isfield(state, 'ks')
	error('kernel structure must be specified (see kernelFactory)');
    end
    if ~isfield(state.ks, 'kernel')
	error('invalid kernel structure');
    end
    if isfield(state, 'learningRate')
	if abs(state.learningRate) > 1
	    error('Learning rate must be less than 1');
	end
    else
	state.learningRate = 0.1; % default learning rate
    end

    state.x = cell(0);
    state.coeff = [];
    state.n = 0;

    out = state;
    return
end

if state.n > 0
    for k = 1:state.n
	kval = state.ks.kernel(state.ks, x, state.x{k}, state.ksize);
    end
    yhat = sum(state.coeff .* kval);
else
    yhat = 0;
end

if nargin < 3 % Prediction mode
    out = yhat;
    return
end

% Training mode
err = y - yhat;
state.x{end+1} = x;
state.coeff(end+1) = state.learningRate * err;
state.n = state.n + 1;
out = state;
