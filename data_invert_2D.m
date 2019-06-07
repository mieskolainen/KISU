%% Invert data 2D loop
%
% call this from data_main.m
%
% mikael.mieskolainen@cern.ch, 2019

close all;

% Inline function [TEST FUNCTION]
func = @(a,b,kappa) a - kappa*b;

for mode = 1:2

% INEL
if      (mode == 1)
    VALS = INEL;

% NSD
elseif (mode == 2)
    VALS = NSD;

% INEL - a*NSD or NSD - a*INEL
elseif (mode == 3 || mode == 4)
    
    if (mode == 3)
    VALS = INEL;
    end
    if (mode == 4)
    VALS = NSD;
    end
    
    % Find matching pair
    for i = 1:length(INEL)
        for j = 1:length(NSD)
            if (INEL{i}.N      == NSD{j}.N       && ...
                INEL{i}.sqrts  == NSD{j}.sqrts   && ...
                INEL{i}.eta(1) == NSD{j}.eta(1) && ...
                INEL{i}.eta(2) == NSD{j}.eta(2))
                
                % Find the minimum number of bins
                D = min(length(INEL{i}.n), length(NSD{j}.n));
                
                alpha_values = linspace(1e-2, 1.0, 1e3);
                
                % Find alpha which keeps the probabilities positive
                alpha = 0.0;
                for k = 1:length(alpha_values)
                    
                    if (mode == 3)
                        ghat = func(INEL{i}.g(1:D), NSD{i}.g(1:D), alpha_values(k)); 
                    end
                    if (mode == 4)
                        ghat = func(NSD{i}.g(1:D), INEL{i}.g(1:D), alpha_values(k));
                    end
                    
                    if (sum(ghat(ghat < 0)) < 0)   % Went negative already
                        alpha = alpha_values(k-1); % Take previous
                        break;
                    end
                end
                % Now subtract
                if (mode == 3)
                    ghat = func(INEL{i}.g(1:D), NSD{i}.g(1:D), alpha); 
                end
                if (mode == 4)
                    ghat = func(NSD{i}.g(1:D), INEL{i}.g(1:D), alpha);
                end
                break; % Found the match loop
            end
        end
        
        VALS{i}.n = VALS{i}.n(1):length(ghat - VALS{i}.n(1)) - 1; % Bins
        VALS{i}.g = ghat / sum(ghat); % Re-normalize
        
        if (mode == 3)
            VALS{i}.type = 'INEL_MINUS_NSD';
            title_str = '$\mathcal{N}$(INEL - $\alpha$NSD)';
        end
        if (mode == 4)
            VALS{i}.type = 'NSD_MINUS_INEL';
            title_str = '$\mathcal{N}$(NSD - $\alpha$INEL)';
        end
    end
end

% Loop over datasets
for k = 1:length(VALS)

type   = VALS{k}.type;
N      = VALS{k}.N;
sqrts  = VALS{k}.sqrts;
eta    = VALS{k}.eta;
n      = VALS{k}.n;
EVENTS = VALS{k}.EVENTS;

g      = VALS{k}.g;
g_pos  = VALS{k}.g_pos;
g_neg  = VALS{k}.g_neg;


%% Find lambda parameter

%param.mu = 3;
%input = g * EVENTS;
%[param.lambda, fig1] = findlambda(input, [], param, lambdas);

%% Inverse

% Technical domain extension factor (2 x enough here, speed improvement)
extension = 2.0;

tic;

% x-axis maximum value
MAX_X = 150;

fig1 = figure('Visible','Off');

mu_val = linspace(0.01, 8.0, 60);
f_hat = zeros(length(g), length(mu_val));

input = g;
input = input * EVENTS;     % Scale with approximate statistics

parfor kk = 1:length(mu_val)
    
    mu = mu_val(kk);
    
    [f_iters, rho, F_tot, g_err, reg_err, g_hat] = ...
        kisusolver(input, [], mu, param.R, param.lambda, param.algo, extension);
    f_hat(:,kk) = f_iters(:,end);                 % Take the last iteration result
    f_hat(:,kk) = f_hat(:,kk) / sum(f_hat(:,kk)); % Normalize to prob density
end

% Note transpose for f_hat
imagesc(n, mu_val, log10(f_hat')); colormap(hot);
set(gca,'YDir','normal');
if (strcmp(type,'INEL') && N == 1)
    colorbar;
end
caxis(log10([10^(-3.5) 10^(-0.5)]));

% INEL or NSD
if (strcmp(type,'INEL') || strcmp(type,'NSD'))
title(sprintf('%s $N_{ch} \\geq %d$: $\\sqrt{s}=%0.1g$ TeV, $\\eta \\in [%0.1f, %0.1f]$', ...
      type, N, sqrts, eta(1), eta(2)), 'interpreter','latex');

% INEL - NSD test
else
    title(sprintf('%s, $\\alpha = %0.2f$, $N_{ch} \\geq %d$: $\\sqrt{s}=%0.1g$ TeV, $\\eta \\in [%0.1f, %0.1f]$', ...
      title_str, alpha, N, sqrts, eta(1), eta(2)), 'interpreter','latex');
end
axis square;

% Tick marks on top of data
set(gca,'Layer','top')

xticks([0:25:MAX_X]);

axis([0 MAX_X min(mu_val) max(mu_val)]);

xlabel('$N_{ch}$',      'interpreter', 'latex');
ylabel('$\mu$', 'interpreter', 'latex');

% Tick marks on top of data
set(gca,'Layer','top')
% ------------------------------------------------------------------------

% Print out
filename = sprintf('2D_%s_N%d_SQRTS_%0.0f_ETAMIN_%0.1f_ETAMAX_%0.1f_SKIP0BIN_%d.pdf', ...
                    type, N, sqrts*1e3, eta(1), eta(2), SKIP0BIN);
print_cmd = sprintf('../lhcfigs2D/%s', filename);
print(fig1, print_cmd, '-dpdf', '-painters');
system(sprintf('pdfcrop --margins 1 ../lhcfigs2D/%s ../lhcfigs2D/%s', filename, filename));

toc;
close all;

%%
end
end

toc;
