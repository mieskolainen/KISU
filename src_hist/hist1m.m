% 1D histogram
%
% input:       X =  Data events (N)
%          xedge =  x-range bin EDGES [xmin, ... , ... , xmax]
%
% output:      h =  Count vector
%
% mikael.mieskolainen@cern.ch, 2017

function h = hist1m(X, xedge)

h = zeros(length(xedge)-1, 1);

% Loop in x-direction
for i = 1:length(xedge)-1
    
    % Collect events
    xmin  = xedge(i);
    xmax  = xedge(i+1);
    
    % Select events within this [xmin xmax] interval
    ind = (xmin <= X & X < xmax);
    
    % Sum events
    h(i) = sum(ind(:));
end

end