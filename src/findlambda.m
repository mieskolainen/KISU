% Loop over lambda regularization parameter values
%
% INPUT:
%                g = Input distribution histogram
%                b = Background histogram, set [] if none
%            param = Algorithm parameter object struct
%          lambdas = Lambda values, e.g. logspace(-2.301, 1.00, 100);
% OUTPUT:
%           lambda = Optimal lambda according to re-projection chi^2 ~ ndf
%             fig1 = Figure handle
%
% Same input as for mainsolver()
%
% mikael.mieskolainen@cern.ch, 2019

function [lambda, fig1] = findlambda(g, b, param, lambdas)

D = length(g);
% Find optimal lambda parameter

reg_errors = zeros(size(lambdas));    % Regularization cost
ghat_errors = zeros(size(lambdas));   % Reconstruction cost
f_hats = zeros(length(g), length(lambdas));

for l = 1:length(lambdas)
    % Calculate estimate
    [f_iters, ~, ~, g_err, reg_err, ~] = ...
        kisusolver(g, b, param.mu, param.R, lambdas(l), param.algo);
    f_hats(:,l) = f_iters(:,end);
    ghat_errors(l) = g_err(end);
    reg_errors(l) = reg_err(end);
end

% Plot lambda curves
fig1 = figure;
a = ghat_errors;% / sum(ghat_errors);
b = reg_errors;% / sum(reg_errors);
c = a + b;

h0 = plot(lambdas, ones(size(lambdas))*D, 'k:', 'color', [ones(3,1)*0.5], 'linewidth', 1.0); hold on; % horizontal line
h1 = plot(lambdas, a, 'r-', 'linewidth', 1.0); 
h2 = plot(lambdas, b, 'k-', 'linewidth', 1.0);
h3 = plot(lambdas, c, 'b--', 'color', [0.4 0.4 0.5], 'linewidth', 1.0); % sqrt(a^2+b^2)

text(lambdas(1)*1.2, D*1.5, '$D$', 'interpreter','latex', 'color', [ones(3,1)*0.5]); % text 'D'

axis tight; axis square;
set(gca, 'xscale','log','yscale','log');
xlabel('$\lambda$', 'interpreter','latex');

%
% I.) Find "discrepancy param.principle" minimum, re-projection chi^2 ~ ndf point
if      (strcmp(param.principle, 'a_minus_ndf_min'))
disc = abs(a - D);
[~,min_ind] = min(disc);
lambda = lambdas(min_ind);

% Plot vertical line indicating minimum
plot(ones(2,1)*lambda, [min([a b])*1e-2 a(min_ind)], 'k-.');

% II.) Find a+b=c minimum, and choose as the lambda
elseif (strcmp(param.principle, 'a_plus_b_min'))

[~,min_ind] = min(c);
lambda = lambdas(min_ind);

% Plot vertical line indicating minimum
plot(ones(2,1)*lambda, [min([a b])*1e-2 c(min_ind)], 'k-.');

% III.) Find minimum according to the a = b crossing point
elseif (strcmp(param.principle, 'a_minus_b_min'))

disc = abs(a-b);
[~,min_ind] = min(disc);
lambda = lambdas(min_ind);

% Plot vertical line indicating minimum
plot(ones(2,1)*lambda, [min([a b])*1e-2 a(min_ind)], 'k-.');

else
   error('findlambda: Unknown param.principle =' + param.principle); 
end

if (lambda > 0.2) % CUTOFF
    lambda = 0.2;
end

% Set legend
leg = legend([h1 h2 h3], {'$\alpha:\, \chi^2_{\hat{g}}$', ...
                          '$\beta:\, \chi^2_{\nabla^2 \hat{f} }$', ...
                          '$\alpha + \beta$'} );
set(leg, 'interpreter','latex', 'location','southeast');
legend('boxoff');

%axis([min(lambdas) max(lambdas) min([a b])/3 max([a b])*3]);

end