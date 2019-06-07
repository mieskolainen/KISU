% Sample N events from a discrete density
% 
% Input:
%         pdf = discrete pdf from which we sample
%           N = Number of events
%   normalize = Normalize the output histogram (true / false)
% 
% Output:
%           y = Output histogram
%
% mikael.mieskolainen@cern.ch, 2019

function y = sampledensity(pdf, N, normalize)

pdf = pdf(:) / sum(pdf);

% Sample with replacement
sample = randsample(1:length(pdf), N, true, pdf);
y = zeros(size(pdf));
for i = 1:length(y)
    y(i) = sum(sample == i);
end
if (normalize) % normalize to 1
   y = y / sum(y); 
end

end