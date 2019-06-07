% Visualize the simulation results
% ------------------------------------------------------------------------
%
% Output from mainsolver()
%
% mikael.mieskolainen@cern.ch, 2019

function figs = visualize(output, g, param, tru, bc)

dX = bc(2) - bc(1); % Binwidth
N = sum(g);         % Number of events
D = length(g);      % NDF

% Calculate the statistics based on the mother sample

% Median
f_hat = prctile(output.F_star', 50)';
f_hat_0 = prctile(output.F_0', 50)';

% Estimate 95CL = 2.5 ... 97.5 confidence levels of our point estimate with bias correction
F_min = prctile(output.F_star', 2.5)';  F_min = F_min + 1e-9;  % Due to log scale
F_max = prctile(output.F_star', 97.5)'; F_max = F_max + 1e-9; 

% Other statistics based on last inner bootstrap

% Error on g (projection)
E_min = prctile(output.E', 16)';
E_max = prctile(output.E', 84)';

% Error on f (estimate)
EF_min = prctile(output.EF', 16)';
EF_max = prctile(output.EF', 84)';

% Regularization cost
RG_min = prctile(output.RG', 16)';
RG_max = prctile(output.RG', 84)';

g_err  = prctile(output.E', 50)';
f_err  = prctile(output.EF', 50)';
RG_med = prctile(output.RG', 50)';

% Get ratio plot
[fig,ax] = ratioplot();
figs{1} = fig;

% Visualization limits
x_limit = bc(end);

% Set the remaining axes properties
set(ax{1},'YScale','log');
axes(ax{1});

% PLOT1 : Main plot
%    xsteps = linspace(0,x_limit,length(f));

plot(bc, tru*N*dX, 'k-.','linewidth', 1); hold on; % True
stephist(bc, g, 'r-');                             % Folded
stephist(bc, f_hat_0, 'g-', 'linewidth',1);        % First iteration
stephist(bc, f_hat, 'b-', 'linewidth',1);          % Unfolded

% Signal to noise ratio (signal times 1.96, because F_max, F_min are 2 sigma
% CLs)
snr = f_hat ./ ((abs(F_max - f_hat) + abs(F_min - f_hat))/2) * 1.96;

% Smooth it
snr = medfilt1(snr,5);
ind = 1;
sigma_cut = 3;
for z = 1:length(snr)
   if (snr(z) < sigma_cut) % Below sigma level
      ind = z; break; 
   end
end
% Sigma cut region
snr_color = ones(3,1)*0.7;
plot(linspace(0, x_limit, 3), ones(3,1)*f_hat(ind), 'color', snr_color, 'linestyle','--');
text(bc(1)*10, f_hat(ind)*1.5, sprintf('$\\hat{f} > %0.0f\\hat{\\sigma}$', sigma_cut), ...
    'interpreter','latex', 'color', snr_color);

% Confidence levels
color = [0 0 1];
stepfill(bc, F_max, F_min, color, color, 0.05); hold on;

set(gca,'yscale','log');
lh1 = legend('True $f$','Convoluted $g$','Estimate $\hat{f}_0$','Bias iterated $\hat{f}_{\star}$'); 
set(lh1,'interpreter','latex','location','northeast');
legend('boxoff');

axis square; 

ylabel(sprintf('counts $ |_{\\Delta x=%0.1f}$', dX),'interpreter','latex');
title(sprintf('$\\mu = %0.2f, N = 10^{%0.0f}, \\lambda = %0.2f, D = %0.0d$', ...
        param.mu, log10(N), param.lambda, D),'interpreter','latex');
axis([0 x_limit 1 10^(log10(N/2))]);

% Remove lowest tick mark, because it can overlap with the lower plot
hold(ax{1},'on');
yTick = get(ax{1},'YTick');
set(ax{1}, 'YTick', yTick(2:end));

% --------------------------------------------------------------------
% PLOT2 : Ratio plot
axes(ax{2});

% horizontal line
plot(linspace(0, x_limit, 10), ones(10,1), 'k-'); hold on;

% Note: No error propagation needed here, because B (= truth) is exact
% thus taylor expansion of ratio f = A/B reduces to what is below
stepfill(bc, (F_max(:) ./dX )./(tru*N), ...
            (F_min(:) ./dX )./(tru*N), color, color, 0.05); hold on;
stepbar(bc, (f_hat_0(:) ./dX )./(tru*N), 'g-'); %
stepbar(bc, (f_hat(:) ./dX )./(tru*N), 'b-'); %

axis([0 x_limit 0.5 1.5]);
xlabel('$x$','interpreter','latex');
ylabel('$\hat{f}/f$','interpreter','latex');

% Set the remaining axes properties
set(ax{2},'XMinorTick','on','YMinorTick','on');


%% Metrics as a function iterations

figs{2} = figure;

plot(0:length(f_err)-1, f_err,'b-'); hold on;
plot(1:length(g_err), g_err,'r-');
plot(1:length(RG_med), RG_med,'k-');

color = [0 0 1]*0.9;
plotfill(0:length(EF_max)-1, EF_max, EF_min, color, color, 0.2); hold on;

color = [1 0 0]*0.9;
plotfill(1:length(E_max), E_max, E_min, color, color, 0.2); hold on;

color = [1 1 1]*0.5;
plotfill(1:length(RG_max), RG_max, RG_min, color, color, 0.2); hold on;

set(gca,'yscale','log','xscale','log');

%plot(f_err,'b-');
xlabel('Iteration $k$','interpreter','latex'); %ylabel('$\chi_k^2$','interpreter','latex');
axis square;
title(sprintf('$\\mu = %0.2f, N = 10^{%0.0f}, \\lambda = %0.2f, D = %0.0d$', param.mu, log10(N), param.lambda, D),'interpreter','latex');
lh = legend(sprintf('$\\chi^2_{\\hat{f}}, \\min = %0.2E$', min(f_err)), ...
            sprintf('$\\chi^2_{\\hat{g}}, \\min = %0.2E$', min(g_err)), ...
            sprintf('$\\chi^2_{\\nabla^2 \\hat{f} }, \\min = %0.2E$', min(RG_med)));
legend('boxoff');
set(lh,'Interpreter','latex','location','northeast')

axis([0 param.R min([f_err(:);RG_med(:);g_err(:)])*0.1 max([f_err(:);RG_med(:);g_err(:)])]);


end
