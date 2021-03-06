% Stepfill with EDGE bins
%
% input:      x = x-axis BIN EDGES (N)
%             u = upper curve data (N-1)
%             l = lower curve data (N-1)
%         color = Fill color
%          edge = Edge color
%  transparency = Fill tranparency [0...1]
%
% mikael.mieskolainen@cern.ch, 2017

function stepfilledge(x,u,l,color,edge,transparency)

if (length(u) ~= length(l))
    error('stepfilledge:: length(u) ~= length(l) !');
end
if (length(u) ~= length(x)-1)
    error('stepfilledge:: x,u,l, input dimensions do not match !');
end

x_ = [];
u_ = [];
l_ = [];

% Create steps
for i = 1:length(x)-1
    delta = (x(i+1) - x(i));

    x_(end+1) = x(i);
    x_(end+1) = x(i) + delta;
    u_(end+1) = u(i);
    u_(end+1) = u(i);
    l_(end+1) = l(i);
    l_(end+1) = l(i);
end

x_(1) = x_(1) + 1e-9; % For logarithmic x-scale, if the left edge is at 0

% Plot each bar separatily, for robustness (sometimes fails otherwise for noisy data)
for i = 1:2:length(x_)
    plotfill(x_(i:i+1),u_(i:i+1),l_(i:i+1),color,edge,transparency); hold on;
end
end