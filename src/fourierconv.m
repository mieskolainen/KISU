% Linear Convolution on Fourier Domain
%
% Faster than conv() for long&long inputs.
% However, conv() faster with short&short, long&short inputs.
%
% mikael.mieskolainen@cern.ch, 2019

function [y,Y] = fourierconv(x, h, N_fix)

N  = length(x)+length(h)-1; % Length of convolution
N2 = pow2(nextpow2(N));     % Take the nearest power of 2 > N, for speed

if (nargin == 3)
    N2 = N_fix;
end

Y = fft(x, N2) .* fft(h, N2); % IFFT( FFT * FFT )
y = real(ifft(Y));
y = y(1:N);                   % Truncate to the original width

end