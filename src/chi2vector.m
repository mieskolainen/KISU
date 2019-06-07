% Calculate Chi^2 for a set of estimates
%
% INPUT:
%       X = D (dim) x N (sample size)
%     tru = Sampled continous pdf (truth)
%
% OUTPUT:
%     err = chi2 value
%
% mikael.mieskolainen@cern.ch, 2019

function err = chi2vector(X, tru)

err = zeros(size(X,2),1);

% Calculate chi2 for all estimates
for i = 1:size(X,2)
    chi2_per_bin = (tru - X(:,i)).^2 ./ tru; 
    chi2_per_bin( tru < 1 ) = 0; % At least one event per bin, null rest
    err(i) = sum(chi2_per_bin);
end

end