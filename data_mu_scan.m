%% Optimal mu-value scan loop
%
% call this from data_main.m
%
% mikael.mieskolainen@cern.ch, 2019

close all;

lambda = 0.0; % Regularization strength
STEPS  = 50; % mu_value steps

% Minimum N value, this is adjusted in a way that we can skip the first
% NBD peak, we are interested in the secondary NBD peak maximum negative
% derivative here
N_min_ind = 10;

sqrts_values = [0.9 7.0 8.0];

for s = 1:length(sqrts_values)
fig1 = figure;
legs = {};

% mu-valus at minimum
min_mu = [];

max_y = -1e99;
min_y =  1e99;

for i = 1:length(NSD)
    if ( NSD{i}.N == 0 && NSD{i}.sqrts == sqrts_values(s) ) %&& NSD{i}.eta(2) == 5.0) 
        
        type   = NSD{i}.type;
        N      = NSD{i}.N;
        sqrts  = NSD{i}.sqrts;
        eta    = NSD{i}.eta;
        n      = NSD{i}.n;
        g      = NSD{i}.g;
        EVENTS = NSD{i}.EVENTS;
        
        % Find value of mu with maximal 
        fprintf('HEP \n');
        
        % Run KISU algorithm
        mu_values = linspace(1.0, 6.0, STEPS);
        dip_values = zeros(length(mu_values),1);
        parfor k = 1:length(mu_values) % ** PARALLEL LOOP HERE **
            
            mu = mu_values(k);
            
            input = g * EVENTS; 
            [f_iters, rho, F_tot, g_err, reg_err, g_hat] = ...
              kisusolver(input, [], mu, param.R, param.lambda, param.algo, extension);
            f_hat = f_iters(:,end); % Take the last iteration
            
            %ind0 = f_hat > 0;
            %entropy = -sum(f_hat(ind0).*log2(f_hat(ind0)));
            
            derivative = diff(f_hat);
            dip_values(k) = min(derivative(N_min_ind:end));
        end
        
        % Spline interpolate to extract sub-sampled minimum
        mu_values_interp1  = mu_values(1):1e-3:mu_values(end);
        dip_values_interp1 = interp1(mu_values, dip_values, mu_values_interp1, 'spline');
        
        [~,min_ind] = min(dip_values_interp1);
        min_mu(end+1) = mu_values_interp1(min_ind);
        
        plot(mu_values, dip_values, 'linewidth', 1.0); hold on;
        xlabel('$\mu$','interpreter','latex');
        ylabel('min $d\hat{f}(\mu)/dN$','interpreter','latex');
        xticks([min(mu_values):0.5:max(mu_values)]);
        axis square;
        legs{end+1} = sprintf('$\\eta \\in [%0.1f, %0.1f]$', eta(1),eta(2));
        title(sprintf('%s $N_{ch} \\geq %d$: $\\sqrt{s} = %0.1g$ TeV', type, n(1), sqrts),'interpreter','latex');
        
        % Save extremum for visualization
        if (min(dip_values_interp1) < min_y), min_y = min(dip_values_interp1); end
        if (max(dip_values_interp1) > max_y), max_y = max(dip_values_interp1); end
    end
end

mu_hat     = mean(min_mu);                       % Use mean
mu_hat_err = std(min_mu) / sqrt(length(min_mu)); % Standard error of mean

x = 3.5;            % coordinates
y = max_y*1.25;
text(x, y, sprintf('$\\hat{\\mu} = %0.2f \\pm %0.2f$', mu_hat, mu_hat_err), 'interpreter', 'latex', 'fontsize', 14);
axis([min(mu_values) max(mu_values) min_y * 1.25 max_y * 0.5]);

l = legend(legs,'location','southeast');
set(l,'interpreter','latex');
legend('boxoff');

filename  = sprintf('mu_scan_%s_N%d_SQRTS_%0.0f_SKIP0BIN_%d.pdf', type, N, sqrts*1e3, SKIP0BIN);
print_cmd = sprintf('../lhcfigs/%s', filename);
print(fig1, print_cmd, '-dpdf', '-painters');
system(sprintf('pdfcrop --margins 1 ../lhcfigs/%s ../lhcfigs/%s', filename, filename));

close all;
end
