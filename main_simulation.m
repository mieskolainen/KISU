% Main simulation of K-fold stochastic autoconvolution inverse
%
% mikael.mieskolainen@cern.ch, 2019

clear; close all;
addpath('./src_hist');
addpath('./src');

system('mkdir ../figs');

PRINT_ON = true;

%%
% ------------------------------------------------------------------------
% Simulation scenarios
N_values     = [1e5 1e6];     % Number of events
mu_values    = [0.5 2.0 5.0]; % Different mu-values
distr_values = [1 2]; %       % Different distributions

% ------------------------------------------------------------------------
% Algorithm parameters

param.algo   = 5;    % Algorithm 5 for Poisson case, 3 for Gaussian approx

param.N_moth = 100;  % Number of mother bootstrap samples (as large as possible)
param.N_boot = 30;   % Number of bootstrap daughter samples (as large as possible)
param.N_bias = 3;    % Number of daughter bootstrap bias iterations (>= 1)
param.R      = 50;   % Number of iterations in the inverse algorithm (>= 50)

% Regularization strength choise principle
param.lambda = 0.1;  % Default regularization parameter

%param.principle = 'a_minus_ndf_min'; % Discrepancy principle
param.principle = 'a_plus_b_min';     % Variational equilibrium
%param.principle = 'a_minus_b_min';   % Test

% Lambda values to scan
lambdas  = logspace(-2.301, 1.00, 30);

% ** Set per simulation below **
param.mu = 0.0;

tic;

%% Loop over test distributions
for distr = distr_values


%% Generate the truth distribution f(X), with sum_x f = 1
if (distr == 1)
    % Negative binomial fit to P(Nch) | eta = [-1...1] at sqrt(s) = 7 TeV
    % https://arxiv.org/pdf/1202.4221.pdf
    
    max_xval = 400; % Needs to be high enough, not to bias the tail
    bins = 1e4;     % This needs to be very high => we are approximating continuous distribution
    samples = linspace(0,max_xval,bins);
    
    % Negative-Binomial distribution
    nbdpdf = @(n,mu,k) gamma(k+n)./(gamma(k).*gamma(n+1)).*(mu./(k+mu)).^n.*(k./(k+mu)).^k;
    f = nbdpdf(samples,12.5,1.4); f(isnan(f)|isinf(f)) = 0;
    f = f(:) / sum(f);
    
    % Binwidth proposal
    dX = 2.0;
    
elseif (distr == 2)
    
    % Exponential pdf
    max_xval = 100; % Needs to be high enough, not to bias the tail
    bins = 1e4;     % This needs to be very high => we are approximating continuous distribution
    samples = linspace(0,max_xval,bins);
    b = 6;          % exp distribution mean
    f = exppdf(samples, b);
    f = f(:) / sum(f);
    
    % Binwidth Proposal
    dX = 2.0;
    
elseif (distr == 3)
    
    % Gamma pdf
    max_xval = 400; % Needs to be high enough, not to bias the tail
    bins = 1e4;     % This needs to be very high => we are approximating continuous distribution
    samples = linspace(0,max_xval,bins);
    f = gampdf(samples, 2, 5);
    f = f(:) / sum(f);
    
    % Binwidth Proposal
    dX = 2.0;
    
else
    error('Unknown simulation distribution %d', distr);
end

% Draw simulated
%stephist(samples,f); axis tight;
%set(gca,'yscale','log');


%% Mapping of a continuous x value in R to a discrete x bin in Z

xdisc = @(xcont) ceil(xcont/max_xval*bins);


%% Now generate MC sample with Poisson pileup

for  N = N_values
for mu = mu_values

DATA = zeros(N,1);

% Create Poisson probabilities of k = 1,2,...,k_max for fast random numbers
% (poissrnd() is very slow)
max_k = 100;
poiss_pdf = poisspdf(1:max_k, mu);
poiss_pdf = poiss_pdf / sum(poiss_pdf);
poissons = randsample(1:max_k, N, true, poiss_pdf);

% Function values
maxf = max(f);

tic;
for i = 1:N
    
    % Draw the pileup number
    k = poissons(i);
    
    % Draw k samples from the distribution f via
    % von Neumann (acceptance-rejectance)
    x_values = zeros(k,1);
    for j = 1:k
        while (true)
            % First draw x-axis value
            xval = rand(1)*max_xval;
            % Draw g-axis value
            if (rand(1)*maxf < f(xdisc(xval)) )
                x_values(j) = xval; break; % Accept
            end
        end
    end
    DATA(i) = sum(x_values);
    
    if (mod(i,1e5) == 0)
       fprintf('Samples %d/%d generated \n', i, N); 
    end
end
toc;

%%
clear poissons;


%% Now discretization of the problem with histogram
% Note that bin width (or number of bins) is the first regularizator of the
% problem, (low-pass filter in spectral domain, basically)

% Compute bins
minimum = 0;
maximum = prctile(DATA, 99.99); % This cuts the upper tail

% By constant linear discretization/binning

% By constant linear bin width
D = round( (maximum-minimum)/dX );
binEdges = linspace(minimum, maximum, D+1);
dX = binEdges(2) - binEdges(1); % Update
%}

%
if (D > 350)
    D = 350;
    binEdges = linspace(minimum, maximum, D+1);
    dX = binEdges(2)-binEdges(1); % Discretization (binwidth)
end
%}

aj = binEdges(1:end-1)';       % Bins lower edge
bj = binEdges(2:end)';         % Bins upper edge
bc = (aj + bj)/2;               % Bins center

% Calculate histogram
g = histcounts(DATA, binEdges)';


%% Create "ground truth" sampled PDF for a comparison
xsteps = linspace(0,max_xval,length(f));

% Find indices corresponding the center points of histogram
f_ind = zeros(size(bc));
for k = 1:length(bc)
    [~,f_ind(k)] = min(abs(xsteps - bc(k)));
end
tru  = f(f_ind) / (xsteps(2)-xsteps(1)); % Normalize with discretization


%% Set mu-parameter

% mu-measurement noise perturbation (for a more realistic simulation)
% mu_measured = mu;  % no noise case, fine for us
sigma_mu = 0.0;      % noise std
mu_measured  = mu + sigma_mu*randn(1)*mu;

param.mu     = mu;   % Compound Poisson parameter


%% Find "optimal" regularization parameter

b = []; % No background in simulation

[param.lambda, fig1] = findlambda(g, b, param, lambdas);

if (PRINT_ON == true)
    filename = sprintf('LAMBDA_D%d_A%d_mu_%0.2f_N_10E%0.0f_%s.pdf', distr, param.algo, mu, log10(N), param.principle);
    % -painters renderer more stable than opengl
    print_cmd = sprintf('../figs/%s', filename); 
    print(fig1, print_cmd, '-dpdf', '-painters');
    system(sprintf('pdfcrop --margins 1 ../figs/%s ../figs/%s', filename, filename));
end
close all;

%% INVERSE SOLUTION

output = mainsolver(g, b, param, tru, bc);


%% Visualize

figs = visualize(output, g, param, tru, bc);

% Inverse distribution
if (PRINT_ON == true)
    filename = sprintf('RESULT_D%d_A%d_mu_%0.2f_N_10E%0.0f_%s.pdf', distr, param.algo, mu, log10(N), param.principle);
    % -painters renderer more stable than opengl
    print_cmd = sprintf('../figs/%s', filename); 
    print(figs{1}, print_cmd, '-dpdf', '-painters');
    system(sprintf('pdfcrop --margins 1 ../figs/%s ../figs/%s', filename, filename));
end

% Error evolution
if (PRINT_ON == true)
    filename = sprintf('ERROR_D%d_A%d_mu_%0.2f_N_10E%0.0f_%s.pdf', distr, param.algo, mu, log10(N), param.principle); 
    % -painters renderer more stable than opengl
    print_cmd = sprintf('../figs/%s', filename);
    print(figs{2}, print_cmd, '-dpdf', '-painters');
    system(sprintf('pdfcrop --margins 1 ../figs/%s ../figs/%s', filename, filename));
end
close all;

end % mu value loop
end % N of events loop
end % Different distributions loop

toc;
