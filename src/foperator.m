% Forward operator F(mu) : f -> g of K-fold compound Poisson autoconvolution
%
% Implementation using FFT
%
% input:    f_hat = input pdf function (X domain)
%           mu    = Poisson mu
%
% output:   g_hat = output pdf function in same domain as input (X domain)
%      g_hat_full = output pdf function in larger domain (Y domain)
%               h = Forward mapping as a convolution kernel vector
%
% mikael.mieskolainen@cern.ch, 2019

function [g_hat,g_hat_full,h] = foperator(f, mu)

% If arbitrary input distribution, not just Poisson mu
if (length(mu) ~= 1) 
   
   [g_hat,g_hat_full] = foperatorToeplitz(f_hat, mu);
   return;
end

% Poisson case =>
events = sum(f);
f = f / sum(f); % Make it probability distribution with sum = 1

% Fourier transform with long tail
meanval = mu/(1-exp(-mu)); % Effective Poisson mu
N = length(f)*ceil((meanval)*5); % 5 x to capture the tail
N = pow2(nextpow2(N));     % Take the nearest power of 2 > N, for speed
F = fft(f, N);

% Map the characteristic function (chf) via Poisson with 0-bin removed
G = 1/(exp(mu)-1)*(exp(mu*F)-1);

% Calculate the corresponding linear kernel
H = G./F; H(isnan(H)) = 0; % in some extremes, division by zero may happen
h = real(ifft(H));
h = h(1:length(f));

% Inverse map with full tail
g_hat_full = real(ifft(G));

% Truncated to the same domain as input
g_hat = g_hat_full(1:length(f))*events;

end