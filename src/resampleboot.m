% Bootstrap multinomial resampling with replacement
%
% input:    x = d x N (d = dimension, N = sample size)
% output:   y = d x N
%
% mikael.mieskolainen@cern.ch, 2019

function y = resampleboot(x)

[~,N] = size(x);

y = x(:, randi(N, [1 N]));

end