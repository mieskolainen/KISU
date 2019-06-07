% Total K-fold Inverse of Stochastic Autoconvolution solver routine
% with bootstrapping for statistical uncertainties
% ------------------------------------------------------------------------
%
% Input:     g   =   The distribution to be inverted
%            b   =   Background histogram, if none put []
%        param   =   Parameter object struct
%          tru   =   Ground truth (for simulations)
%           bc   =   Histogram bincenters (constant histogram binning)
%
% mikael.mieskolainen@cern.ch, 2019

function output = mainsolver(g, b, param, tru, bc)

dX = bc(2)-bc(1); % Binwidth fixed
N = sum(g);       % Number of events
D = length(g);    % Number bins

% The statistics to be plotted
F_star  = zeros(length(g), param.N_moth);
F_0     = zeros(length(g), param.N_moth);

% Mother Bootstrap
for q = 1:param.N_moth
    
    % Inverse algorithm:
    % Outer loop: Recursive Bias correction
    % Inner loop: Bootstrapping
    
    % MC SAMPLING: Pick Poisson sample with mean given by g bin by bin
    g_star = poissrnd(g);
    
    % Run KISU-solver
    [f_iters, ~, ~, ~, ~, g_hat_i] = ...
        kisusolver(g_star, b, param.mu, param.R, param.lambda, param.algo);
    f_hat_0 = f_iters(:,end);
    f_hat_i = f_hat_0;
    
    % Daughter Bootstrap - Bias correction
    for z = 1:param.N_bias
        
        F  = zeros(length(g), param.N_boot);  % Solutions per bootstrap
        
        G  = zeros(length(g), param.N_boot);  % Re-projection per bootstrap
        E  = zeros(param.R, param.N_boot);    % Errors per bootstrap
        EF = zeros(param.R+1, param.N_boot);  % +1 because of iteration 0
        RG = zeros(param.R, param.N_boot);    % Regularization cost
        
        % Bootstrap sample loop (PARALLEL LOOP)
        parfor i = 1:param.N_boot
            fprintf('Mother sample %d/%d : Daughter bootstrap sample %d/%d [bias iter = %d/%d] \n', ...
                q, param.N_moth, i, param.N_boot, z, param.N_bias);
            
            % MC SAMPLING: Resample data using multinomial sampling
            %g_this  = histcounts(resampleboot(DATA(:)'), binEdges);
            
            % MC SAMPLING: Resample data using fast Poisson sampling
            g_this = poissrnd(g_hat_i);
            
            % RUN KISU-ALGORITHM
            [f_iters, ~, ~, g_err, reg_err, g_hat] = ...
               kisusolver(g_this, b, param.mu, param.R, param.lambda, param.algo);
            F(:,i) = f_iters(:,end);
            G(:,i) = g_hat;
            E(:,i) = g_err;
            RG(:,i) = reg_err;
            
            % 0-th iteration
            X = zeros(D,size(f_iters,2)+1);
            X(:,1) = g_this;
            % 1...R iterations
            for j = 1:size(f_iters,2)
               X(:,j+1) = f_iters(:,j); 
            end

            % Calculate the finite support, in order not to compare
            % outside the sane domain
            support = (tru*dX*N > (1/N) );
            EF(:,i) = chi2vector(X(support,:), tru(support)*N*dX);
        end
        
        % Calculate (iterative) bootstrap bias correction
        
        % By definition: Bias[\hat{theta}] = E[\hat{theta}] - theta
        % Now replace E[.] with bootstrap mean
        bias_f_hat = mean(F,2) - f_hat_i;
        
        % Bootstrap bias corrected point estimate estimate with
        % non-negativity enforced
        f_hat_i = f_hat_0 - bias_f_hat;
        f_hat_i = f_hat_i .* (f_hat_i > 0);
        
        % Calculate new g_hat_i by forward operator. This will be used to
        % generate next bias iteration bootstrap sample
        g_hat_i =  foperator(f_hat_i, param.mu);

    end
    
    % Final point estimate of this mother iteration
    F_star(:,q) = f_hat_i;
    F_0(:,q) = f_hat_0;
end

output.F_star = F_star;    % set of bias iterated estimates
output.F_0 = F_0;          % set o non-bias iterated estimates

% Take last values of these
output.EF = EF;            % true comparison (only simulation)
output.RG = RG;            % regularization cost
output.E = E;              % errors

end