% Stephist with CENTER bins
%
% INPUT: x = Bin center  (N)
%        n = Bin content (N)
%
% One can use together with hist(data, bincenters)
%
% mikael.mieskolainen@cern.ch, 2017

function h = stephist(x, n, varargin)

x = x(:)';
n = n(:)';

if (length(n) ~= length(x))
    error('stephist:: x,n input dimensions do not match !');
end

x(1) = x(1) + 1e-9; % For logarithmic x-scale, if the left edge is at 0
h = stairs([x(1)-(x(2)-x(1))/2 ...
            x-(x(2)-x(1))/2 ...
            x(length(x))+(x(2)-x(1))/2], [0 n 0], varargin{:});
end