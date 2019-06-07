% KISU: K-fold Inverse of Stochastic Autoconvolution
% ------------------------------------------------------------------------
% Note that no estimation of mu is possible naively, because when mu->0
% chi2 (of reprojection) -> 0. Thus, mu must be measured or set by theory.
%
% INPUT:
%          g = Measured histogram   (event counts vector)
%          b = Background function, (event counts vector) or [] for empty
%         mu = Poisson mu, compound parameter (scalar)
%              or a vector of probabilities with sum = 1 (custom distribution)
%          R = Number of iterations     (order of 100 ... 300)
%     lambda = Regularization parameter (scalar, e.g. 0.01)
%       algo = [1,2], 3, [4], 5, [6,7]
%              Use algorithm 3 (gaussian noise) or 5 (poisson counting)
%              the rest are just for debugging and tests. 
%        ext = Domain extension factor (leave out for automatic)
%
% OUTPUT:
%    f_iters = inverse results for each iteration    (matrix of vectors)
%        rho = zero-suppressed compound process mean (scalar)
%      F_tot = autoconvolution operator  (matrix)
%      g_err = re-projection chi2        (scalar)
%    reg_err = regularization chi2       (scalar)
%      g_hat = re-projected distribution (vector)
%
% Gradient of ||del u||^2 is 
% del^Tdel u = div(del u) = -Laplacian u
% 
% https://en.wikipedia.org/wiki/Wolfe_conditions#Armijo_rule_and_curvature
%
% mikael.mieskolainen@cern.ch, 2019

function [f_iters, rho, F_tot, g_err, reg_err, g_hat] = kisusolver(g, b, mu, R, lambda, algo, ext)

% Substract background
if (length(g) == length(b))
    g = g - b;
end

% Zero handling, because the update algorithm is multiplicative, it cannot
% update zeros to non-zeros, thus we add tiny constant
EPSILON = 1e-9;
g(g < EPSILON) = EPSILON;

D = length(g);
f_iters = zeros(D,R);

p = [];
% Poisson case, we use the characteristic function map
if (length(mu) == 1)
    rho = mu/(1-exp(-mu)); % zero suppressed mean
    K = ceil(rho); % Convolution domain extension factor
     
% Arbitrary weighting distribution
else
    K = length(mu); % Convolution domain extension factor
    p = mu(:); % Custom probability distribution for pileup 1,2,3,...
    rho = sum(p(:)' .* (1:length(p)) );
    
    % Operator matrices 
    F_tot = [];
    if (algo >= 2)
        F_tot = eye(D*K);
        F = zeros(K, D*K, D);
        for n = 1:K
            F(n,:,:) = eye(D*K, D);
        end
    end
    % Convoluted vectors here
    c = zeros(D*K, K);
end

% Fixed domain extension factor
if (nargin == 7)
    K = ext; % x times
end
fprintf('Mean PU = %0.2f, Domain extension factor K = %0.1f, Regularization lambda = %0.2E \n', rho, K, lambda);

% Measurement (domain extended K x times)
g_ = zeros(round(D*K),1);
g_(1:length(g)) = g;

% Input as the iter k = 0 initial estimate
f_hat = g_;

% Reconstruction ||g - g_hat|| and regularization ||Reg(x)|| error
g_err   = zeros(R,0);
reg_err = zeros(R,0);

% Regularization (finite difference) matrices
% Note the binwidth 1/dX
GRAD = ngrad(length(f_hat), 1); % Default boundary conditions
DD = GRAD'*GRAD;                % LAPLACIAN (second derivative matrix)

% Modify boundary conditions
DD(1,:) = 0; % Make zeros for first and last bins
DD(end,:) = 0;

% Iter k = 1,2,3...
for k = 1:R
    
    % Forward projection step
    % --------------------------------------------------------------------
            
    % Poisson case, characteristic function transform in Fourier domain
    if (isempty(p))
        [g_hat, ~, h] = foperator(f_hat, mu);
        F = convmatrix(h,length(f_hat));
        F_tot  = F(1:length(f_hat),:);
        g_hat = F_tot*f_hat;

    else
    % Arbitrary case, a sequence of convolution matrices
        c(:,1) = f_hat / sum(f_hat); % NOTE THE NORMALIZATION HERE
        for n = 2:K
            F(n,1:D*n-1,1:D) = convmatrix(c(1:D*(n-1),n-1), D); % SLOW matrix representation
            %c(1:D*n-1,n) = fourierconv(c(1:D*(n-1),n-1), c(1:D,1));
            c(1:D*n-1,n) = conv(c(1:D*(n-1),n-1), c(1:D,1));

            %fprintf('Product : \n'); % SLOWER THAN USING CONV
            %tic;
            %c(1:D*K,n) = squeeze(F(n, :,:))*f_hat;
            %toc;
        end
        % Folded estimate: (Same as the for-loop below, but faster 3x)
        F_box = squeeze(sum(bsxfun(@times, F(:,:,1:D), p), 1));
        F_tot(:,1:D) = F_box;
        g_hat = F_tot*f_hat;

        %}
        % Folded estimate:
        %F_tot = sparse(zeros(D*K));
        %for n = 1:K
        %    F_tot = F_tot + squeeze(F(n,:,:))*p(n);
        %end
    end

    % Inversion step (iterative reconstruction)
    % --------------------------------------------------------------------

    % Algorithm (simple multiplicative Poisson update)
    if (algo == 0)

        % Update rule
        f_hat = (g_ ./ max(g_hat, 1e-12)) .* f_hat;
        
        % Regularization
        f_hat = f_hat - lambda*(DD*f_hat);
        
        % To constrain non-negative results (negative might result from
        % regularization step)
        f_hat(f_hat < EPSILON) = EPSILON;
        
    % Algorithm (simple subtractive Gaussian)
    elseif (algo == 1)
        
        omega = 0.5; % Relaxation parameter
        f_hat = f_hat + omega*( g_ - g_hat);
        
        % Regularization
        f_hat = f_hat - lambda*(DD*f_hat);
        
        % To constrain non-negative results (negative might result from
        % regularization step)
        f_hat(f_hat < EPSILON) = EPSILON;
       
    % Gaussian Kaczmarz ART ("Algebraic Reconstruction Technique")
    elseif (algo == 2)
        
        omega = 0.5; % Relaxation parameter
        
        % Kaczmarz scan of the rows
        for i = 1:size(F_tot,1) 
            ind  = mod(size(F_tot,1), k) + 1; % Index
            ai = F_tot(ind,:);                % Row vector
            s  = omega*(g_(ind) - ai*f_hat)/norm(ai)*ai';
            f_hat = f_hat + s;
            
            % Regularization
            f_hat = f_hat - lambda*(DD*f_hat);
            
            % Non-negative constraint projection
            f_hat(f_hat < EPSILON) = EPSILON;
        end
    
    % Algorithm subtractive Gaussian iteration with Matrix Update
    % https://en.wikipedia.org/wiki/Landweber_iteration
    elseif (algo == 3)
        
        omega = 0.5; % Relaxation parameter
        
        % Landweber type update
        f_hat = f_hat - omega*F_tot'*(g_hat - g_);
        
        % Regularization
        f_hat = f_hat - lambda*(DD*f_hat);
        
        % Non-negative constraint projection
        f_hat(f_hat < EPSILON) = EPSILON;
        
    % Algorithm - Gaussian Conjugate Gradient
    % https://en.wikipedia.org/wiki/Conjugate_gradient_method
    elseif (algo == 4)
        
        % Repeat
        if (z == 1)
            r  = g_ - g_hat; 
            pp = r;
        end

        % Dot 1
        alpha = dot(r,r) / dot(pp, F_tot*pp);

        % Update solution
        f_hat = f_hat + alpha*pp;
        f_hat(f_hat < 0) = 0;

        % Update r and beta
        r_next = r - alpha*F_tot*pp;
        beta = dot(r_next,r_next) / dot(r,r);

        % Update pp
        pp = r +  beta*pp;
        r = r_next;
        
        % Regularization
        f_hat = f_hat - lambda*(DD*f_hat);
        
        % To constrain non-negative results (negative might result from
        % regularization step)
        f_hat(f_hat < EPSILON) = EPSILON;
    
    % Algorithm (Richardson-Lucy/EM with regularization)
    % Non-negativity constraint is by construction
    elseif (algo == 5)
        
        % Update with fidelity
        f_hat = (f_hat ./ max(F_tot'*ones(size(f_hat)), 1e-12)) .* (F_tot'*(g_ ./ max(g_hat, 1e-12) ));

        % Regularization as a separate step
        f_hat = f_hat - lambda*DD*f_hat;
        
        % To constrain non-negative results (negative might result from
        % regularization step)
        f_hat(f_hat < EPSILON) = EPSILON;
        
    % Combined variational gradient
    elseif (algo == 6)
        
        % Proximal operators
        % http://www.seas.ucla.edu/~vandenbe/236C/lectures/proxop.pdf
        
        % Fidelity gradient
        grad1 = F_tot'*((g_hat - g_) ./ max(g_hat, 1e-12));
        
        % Regularity gradient
        grad2 = -lambda*DD*f_hat;
        
        % Update
        gamma = 1.0;
        f_hat = f_hat - gamma*(grad1 + grad2);
        
        % To constrain non-negative results (negative might result from
        % regularization step)
        f_hat(f_hat < EPSILON) = EPSILON;
        
    % Algorithm Direct Tikhonov (SLOW)
    elseif (algo == 7)
        
        % Update with fidelity + lambda * regularization
        f_hat = (F_tot'*F_tot + lambda*DD) \ F_tot*g_;
        
        % To constrain non-negative results (negative might result from
        % regularization step)
        f_hat(f_hat < EPSILON) = EPSILON;
    else
       error(sprintf('itersolver:: Unknown algorithm: %d', algo));
    end
    
    % Calculate costs
    
    % 1. Re-projection chi^2
    g_err(k) = calchi2(g_hat(1:D), g_(1:D));
    
    % 2. Regularization chi^2
    vec = (DD*f_hat ./ sqrt(f_hat));
    vec(f_hat < 1) = 0;
    reg_err(k) = norm(vec)^2;
    
    if (mod(k,50) == 0 || k == 1)
    fprintf('Algo = %d, k = %04d : sum(f_hat) = %0.3f, Re-proj chi2 = %0.1E, Reg-err = %0.1E, NDF = %d \n',  ...
        algo, k, sum(f_hat), g_err(k), reg_err(k), length(g) );
    end
    
    % Output of this iteration
    f_iters(:,k) = f_hat(1:D);
end

% Re-projected
g_hat = g_hat(1:D);

% Make it non-negative
f_iters(f_iters < 0) = 0;

end

% Function to calculate chi2 for counts
function chi2 = calchi2(est, orig)

    chi2_per_i = (est - orig).^2 ./ orig;
    chi2_per_i( orig < 1 ) = 0; % at least one event per bin, null rest
    chi2 = sum(chi2_per_i);
end
