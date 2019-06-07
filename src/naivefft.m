% Naive FFT based K-fold autoconvolution inverse
%
% INPUT:       
%              g = Autoconvoluted histogram values
%         x_vals = x-values vector
%              K = Fixed folding number (1,2,3,...)
%         CUTOFF = Spectral cutoff [0 ... 1]
%
% OUTPUT:  
%          f_hat = Inverse
%
% mikael.mieskolainen@cern.ch, 2019

function f_hat = naivefft(g, x_vals, K, CUTOFF)

N = sum(g);      % Total event count
g_ = g / sum(g); % Normalize to one
G  = fft(g_);

% Take K-th complex root (this will fail as a general strategy)
branch = 0;
FF = complexroot(G,K,branch);
%FF = complexroot(G,mu/(1-exp(-mu)),0);
%FF = complexlog(exp(mu-1)*G + 1);

% Shift to symmetric
f_ = fftshift(FF);

%  Box spectral filter high frequency cutoff [0 ... 1]
L = round(length(f_)*CUTOFF);
f_(1:L) = 0;
f_(end-L+1:end) = 0;

% Shift back
FF = ifftshift(f_);

% Inverse FFT solution
f_hat = real(ifft(FF));
f_hat = f_hat(1:length(x_vals)) * N; % Scale with event count

% Truncate to positive
f_hat(f_hat < 0) = 1e-15;

end