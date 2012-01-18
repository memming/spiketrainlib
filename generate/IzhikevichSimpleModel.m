function st = IzhikevichSimpleModel(currentInjection, neuronType)
% Simulate Eugene M. Izhikevich's simple neuron model to obtain spike trains.
% Input
%   currentInjection: injected current sampled at ms
%   neuronType: 1~20 (7 is class 1, 8 is class 2)
% 
% Originally written by Izhikevich, modified by Memming

pars = [0.02      0.2     -65      6       14 ;... % tonic spiking
        0.02      0.25    -65      6       .5 ;... % phasic spiking
        0.02      0.2     -50      2       15 ;... % tonic bursting
        0.02      0.25    -55     0.05     .6 ;... % phasic bursting
        0.02      0.2     -55     4        10 ;... % mixed mode
        0.01      0.2     -65     8        30 ;... % spike frequency adaptation
        0.02      -0.1    -55     6        0  ;... % Class 1
        0.2       0.26    -65     0        0  ;... % Class 2
        0.02      0.2     -65     6        7  ;... % spike latency
        0.05      0.26    -60     0        0  ;... % subthreshold oscillations
        0.1       0.26    -60     -1       0  ;... % resonator
        0.02      -0.1    -55     6        0  ;... % integrator
        0.03      0.25    -60     4        0  ;... % rebound spike
        0.03      0.25    -52     0        0  ;... % rebound burst
        0.03      0.25    -60     4        0  ;... % threshold variability
        1         1.5     -60     0      -65  ;... % bistability
          1       0.2     -60     -21      0  ;... % DAP
        0.02      1       -55     4        0  ;... % accomodation
       -0.02      -1      -60     8        80 ;... % inhibition-induced spiking
       -0.026     -1      -45     0        80];    % inhibition-induced bursting

a = pars(neuronType, 1);
b = pars(neuronType, 2);
c = pars(neuronType, 3);
d = pars(neuronType, 4);
% I = pars(neuronType, 5);

nT = length(currentInjection);
v = -65;
u = b * v;                 % Initial values of u
st = [];             % spike timings
for t = 1:nT
  if v >= 30
      v = c;
      u = u + d;
      st = [st t];
    end
I = currentInjection(t);
  v = v + 0.5*(0.04*v.^2 + 5*v + 140 -u + I); % step .5 ms (numerical stability)
  v = v + 0.5*(0.04*v.^2 + 5*v + 140 -u + I); 
  u = u + a.*(b.*v-u);
end

st = st / 1000; % convert time to ms
