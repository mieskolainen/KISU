% 2D histogram
%
% input:       X =  Data events (Nx2)
%          xedge =  x-range bin EDGES [xmin, ... , ... , xmax]
%          yedge =  y-range bin EDGES [ymin, ... , ... , ymax]
%
% output:      H =  Count matrix, (1,1) is the origin in (x,y) plane
%
% mikael.mieskolainen@cern.ch, 2017

function H = hist2m(X, xedge, yedge)

H = zeros(length(xedge)-1, length(yedge)-1);

% Loop in x-direction
for i = 1:length(xedge)-1
    
    % Collect events
    xmin  = xedge(i);
    xmax  = xedge(i+1);
    
    % Loop over y-direction
    for j = 1:length(yedge)-1
        
        ymin = yedge(j);
        ymax = yedge(j+1);
        
        % Select events within this [xmin xmax] x [ymin ymax] square
        ind = (xmin < X(:,1) & X(:,1) <= xmax) & (ymin < X(:,2) & X(:,2) <= ymax);
        
        % Sum events
        H(i,j) = sum(ind(:));
    end
end

end