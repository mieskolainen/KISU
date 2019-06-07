% Forward operator F(mu) : f -> g of K-fold autoconvolution
%
% Implementation using TOEPLITZ matrices for arbitrary input distributions.
%
% INPUT:    
%           f_hat = input function (X domain)
%           mu    = Poisson mu (OR arbitrary input distribution)
% OUTPUT:
%           g_hat = output function in same domain as input (X domain)
%      g_hat_full = output function in larger domain (Y domain)
%
% mikael.mieskolainen@cern.ch, 2019

function [g_hat,g_hat_full] = foperatorToeplitz(f_hat, mu)

f_hat = f_hat(:);
D = length(f_hat);

% Poisson case
if (length(mu) == 1)
    
    ethresh = 1e-7; % FIXED
    
    % Construct Poisson probabilities
    p = poisspdf(1:300,mu);
    p = p / sum(p);
    % Find based on threshold, how many probability terms we want to consider
    % in the series
    for i = length(p):-1:1
       if p(i) > ethresh
          K = i; break; 
       end
    end
    p = poisspdf(1:K,mu);
    p = p(:) / sum(p);
    
    % Custom probability distribution for pileup 1,2,3,...
else
    K = length(mu);
    p = mu(:); 
end

% Mean pileup, "zero supressed mean", = mu/(1-exp(-mu)) for Poisson
rho = sum(p(:)' .* (1:length(p)) );

F_tot = eye(D*K);
FF = zeros(K, D*K, D);
for n = 1:K
    FF(n,:,:) = eye(D*K, D);
end

c = zeros(D*K, K);
c(1:D,1) = f_hat / sum(f_hat); % NOTE THE NORMALIZATION HERE
for n = 2:K
    FF(n,1:D*(n)-1,1:D) = convmatrix(c(1:D*(n-1),n-1), D); % SLOW matrix representation
    c(1:D*n-1,n) = conv(c(1:D*(n-1),n-1), c(1:D,1));
end

% Same as above for-loop, but faster 3x
F_box = squeeze(sum(bsxfun(@times, FF(:,:,1:D), p), 1));
F_tot(:,1:D) = F_box;

padded = [f_hat; zeros(D*K-D,1)]; % Padd zeros
g_hat_full = F_tot*padded;
g_hat = g_hat_full(1:D);

end