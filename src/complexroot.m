% n-th Complex root of variable z, k-th branch
%
% mikael.mieskolainen@cern.ch, 2019

function y = complexroot(z,n,k)
    z = z(:);
    n = n(:);
    k = k(:);
    phi = atan2(imag(z),real(z));
    y = abs(z).^(1./n) .* exp(1i*(phi + 2*pi*k)./n);
end