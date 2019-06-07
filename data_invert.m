%% Invert data
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

if (N_bootstrap < 10)
    error('N_bootstrap < 10!');
end

% x-axis maximum value
MAX_X = 150;

legs = {};

% Solve the inverse for each mu-hypothesis
p = [];

% Get colors (7 unique in default map)
linecolors = lines(7);
linecolors(3:4,:) = [linecolors(4,:); linecolors(3,:)]; % Switch order of two

% Create plots
[fig1, ax] = ratioplot();

% ------------------------------------------------------------------------
% Tag to lower plot
axes(ax{2});
%
% Plot horizontal line
ph = plot(linspace(0,MAX_X,2), ones(2,1), 'k-'); hold on;

% Plot data errors
transparency = 0.1;
stepfill(n, (g+g_pos)./g, (g-g_neg)./g, [0 0 0], [1 1 1], transparency);
% ------------------------------------------------------------------------

ii = 1;
for mu = mu_vals
    
    % Bootstrap loop
    f_hat = zeros(length(g), N_bootstrap);
    
    % Mean values
    g_hat_0 = [];

    a = (g - g_neg);
    b = (g + g_pos);
    
    for kk = 1:size(f_hat,2)
         
        if (kk == 1)            % Central value
            input = g * EVENTS; % Scale with approximate statistics
        
        elseif (mod(kk,2) == 0) % Negative boundary
            input = poissrnd(a * EVENTS);
            
        else                    % Positive boundary
            input = poissrnd(b * EVENTS);
            
            % Uniform boundaries turned into Gaussian equivalent
            %{
            sigma = sqrt( (b-a).^2 / 12 ); % From uniform distribution
            input = g + sigma .* randn(length(g), 1); % Gaussian
            %}
        end
        
        fprintf('<Boostrap %d/%d> \n', kk, size(f_hat,2));
        
        [f_iters, rho, F_tot, g_err, reg_err, g_hat] = ...
            kisusolver(input, [], mu, param.R, param.lambda, param.algo, extension);
        f_hat(:,kk) = f_iters(:,end);                 % Take the last iteration result
        f_hat(:,kk) = f_hat(:,kk) / sum(f_hat(:,kk)); % Normalize to prob density
        
        % Re-projection
        if (kk == 1)
            g_hat_0 = g_hat;
        end
    end
    
    axes(ax{1});

    % Plot systematic contours (estimate 95CL = 2.5 ... 97.5 confidence levels)
    stepfill(n, prctile(f_hat', 97.5), prctile(f_hat', 2.5), linecolors(ii,:), [1 1 1], 0.15); hold on; % error contour (no legend for this!)
    
    % Central value
    p(end+1)    = stepbar(n, f_hat(:,1), 'linewidth', 1.0, 'color', linecolors(ii,:));
    
    % Create legend
    legs{end+1} = sprintf('$\\hat{f}, \\mu = %0.1f$', mu);
    
    % --------------------------------------------------------------------
    % Tick to lower plot
    axes(ax{2});
    
    % Plot re-projections / measured
    ratio = (g_hat_0/sum(g_hat_0)) ./ g;
    stepbar(n, ratio, 'linewidth', 1.0, 'color', linecolors(ii,:));
    % --------------------------------------------------------------------
    ii = ii + 1;
end

% Plot data
axes(ax{1});

% error contour (no legend for this!)
stepfill(n, g+g_pos, g-g_neg, [0 0 0], [1 1 1], transparency);
p(end+1) = stepbar(n,g,'k-','linewidth', 1.0); hold on;
legs{end+1} = '$g$ (ALICE)';
set(gca,'yscale','log');
axis([0 MAX_X 1e-4 0.1]);

yticks([1e-3 1e-2 1e-1]);
yticklabels({'10^{-3}', '10^{-2}', '10^{-1}'});

ylabel('$P(N_{ch})$', 'interpreter', 'latex');
l = legend(p, legs);
set(l,'interpreter','latex'); legend('boxoff');

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

% ------------------------------------------------------------------------
% Tick to lower ratio plot
axes(ax{2});

xticks([0:25:MAX_X]);
axis([0 MAX_X 0.5 1.5]);

xlabel('$N_{ch}$',      'interpreter', 'latex');
ylabel('$\hat{g} / g$', 'interpreter', 'latex');

% Tick marks on top of data
set(gca,'Layer','top')
% ------------------------------------------------------------------------

% Print out
filename = sprintf('%s_N%d_SQRTS_%0.0f_ETAMIN_%0.1f_ETAMAX_%0.1f_SKIP0BIN_%d.pdf', ...
                    type, N, sqrts*1e3, eta(1), eta(2), SKIP0BIN);
print_cmd = sprintf('../lhcfigs/%s', filename);
print(fig1, print_cmd, '-dpdf', '-painters');
system(sprintf('pdfcrop --margins 1 ../lhcfigs/%s ../lhcfigs/%s', filename, filename));

toc;
close all;

%%
end
end

toc;

