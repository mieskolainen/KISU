% "Toeplitz matrix or diagonal-constant matrix"
% <https://en.wikipedia.org/wiki/Toeplitz_matrix>
% 
% Input as MATLAB toeplitz function
%
% mikael.mieskolainen@cern.ch, 2019

function T = toeplitzmat(c,r)

% Input check
if (nargin == 1), r = c; end

[c_rows, c_cols] = size(c);
[r_rows, r_cols] = size(r);
if ((c_rows ~= 1 && c_cols ~= 1) || (r_rows ~= 1 && r_cols ~= 1))
    error('toeplitzmat: input should be vectors');
end

% Normal tranpose if necessary
if (c_cols ~= 1), c = c.'; end
if (r_rows ~= 1), r = r.'; end
if (r(1) ~= c(1))
    warning('toeplitzmat: diagonal conflict');
end

% Single (complex) vector input
if (nargin == 1)
    c    = conj(c);
    c(1) = conj(c(1)); % first element back
end

% Construct Toeplitz matrix
N_c = length(r); % note r
N_r = length(c); % note c
T   = zeros(N_r, N_c);

for i = 1:min(N_c, N_r)
    T(i:N_r, i)   = c(1:N_r-i+1);
end
for i = 1:min(N_r, N_c - 1)
    T(i, i+1:N_c) = r(2:N_c-i+1);
end

end