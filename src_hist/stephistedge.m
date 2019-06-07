% Stephist with EDGE bins
% 
% input:   x  =  Bin edges    (N)
%          n  =  Bin content  (N-1)
%     varagin =  Extra arguments 
%
% mikael.mieskolainen@cern.ch, 2017

function h = stephistedge(x, n, varargin)
x = x(:)';
n = n(:)'; n = [n 0];

if (length(n) ~= length(x))
    error('stepfilledge:: x,n input dimensions do not match !');
end

x(1) = x(1) + 1e-9; % For logarithmic x-scale, if the left edge is at 0
h = stairs(x, n, varargin{:});

end
