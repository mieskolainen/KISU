% Filled plots
%
% INPUT:      x = Data points on x-axis
%            up = Upper curve
%           low = Lower curve
%         color = Fill color
%          edge = Edge color
%  transparency = Fill tranparency [0...1]
%
% mikael.mieskolainen@cern.ch, 2017

function h = plotfill(x, up, low, color, edge, transparency)

% Make it row vector
x   = x(:)';
up  = up(:)';
low = low(:)';

% These too
color = color(:)';
edge  = edge(:)';

if (length(up) ~= length(low) || length(low) ~= length(x))
    error('plotfill:: Error, vectors with different lengths!');
end

h = fill([x, fliplr(x)], [up, fliplr(low)], color);
set(h, 'EdgeColor', edge, 'FaceAlpha', transparency, 'EdgeAlpha', transparency);

end