% Discretized gradient operator (d-1) x d operating on vectors
% ------------------------------------------------------------------------
%
% type = 1 (Dirichlet boundary conditions)
%        2 (von Neumann boundary conditions)
%
% Discretized divergence operator is the adjoint of grad
% as -D', such that <grad f, g> = -<f, div g>
%
% mikael.mieskolainen@cern.ch, 2019

function D = ngrad(n, type)

% Create difference matrix, such as for n = 3 and type = 1:
% [-1 1 0
% [0 -1 1]
D = -eye(n);
for i = 1:n-1
    D(i,i+1) = 1;
end
%D = D(1:end-1,:);
D(end,end) = 1;

% if type == 1 (default, Dirichlet boundary conditions)

if (type == 2) % von Neumann type boundary conditions
    D(end,end) = 0;
end

end