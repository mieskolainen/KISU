% Convolution matrix representation of convolution
% as convmat() function
% 
% INPUT:
%         a  =  Vector a to be convoluted
%         n  =  Length of convolution
% OUTPUT:
%         b  =  Convolution matrix
%
% mikael.mieskolainen@cern.ch, 2019

function b = convmatrix(a,n)

[rows, cols] = size(a);
if ((rows ~= 1) && (cols ~= 1))
    error('convmatrix: input a should be vector');
end

% Get using Toeplitz matrix
b = toeplitzmat([a(:); zeros(n-1,1)], [a(1); zeros(n-1,1)]);
if (cols > rows), b = b.'; end % Normal transpose

end