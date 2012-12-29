clear

T = 2;
options.ks = kernelFactory('schoenberg', T, 'gaussian');
options.ksize = [0.1 1];
%options.ksize = options.ks.autoParam(options.ks, sts);
options.learningRate = 0.5;

state = stklms(options);
x = [];
y = 0.3;
state = stklms(state, x, y)
yhat = stklms(state, x)
