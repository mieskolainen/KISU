% Plot K-fold autoconvolution for some analytic functions
%
% mikael.mieskolainen@cern.ch, 2019

clear; close all;
addpath('./src_hist');
addpath('./src');

system('mkdir ../paperfigs');


%% Gamma pdf + Gaussian pdf

%{
max_xval = 400; % Needs to be high enough, not to bias the tail
bins = 5e2;     % This needs to be very high => we are approximating continuous distribution
samples = linspace(0,max_xval,bins);
f = gampdf(samples, 2, 5) + 2.0*normpdf(samples, 40, 15);
f = f(:) / sum(f);
%}


%% Negative Binomial Distribution

max_xval = 400; % Needs to be high enough, not to bias the tail
bins = 1e3;     % This needs to be very high => we are approximating continuous distribution
samples = linspace(0,max_xval,bins);

mynbdpdf = @(n,mu,k) gamma(k+n)./(gamma(k).*gamma(n+1)).*(mu./(k+mu)).^n.*(k./(k+mu)).^k;
f = mynbdpdf(samples,12.5,1.4); f(isnan(f)|isinf(f)) = 0;
f = f(:) / sum(f);
%}


%% Plot true first
f2h = figure;

figure(f2h);
delta = samples(2) - samples(1); % discretization width
semilogy(samples, f / delta, 'k--'); hold on;


%% PLot K-fold convoluted with varying mu-values
%
% REMEMBER to check the integrals!

mu = [0.02 0.1 1.0 2.0 5.0 10.0];

for i = 1:length(mu)
    %N = 1e6; % Controls the precision in the forward operator
    %[g,g_full] = foperator(f*N, mu(i));
    %fprintf('mu = %0.3f , sum(g) = %0.3f, integral = %0.3f \n', mu(i), sum(g), sum(g) / N);
    %g = g / N / (samples(2) - samples(1)); % Make it a proper continuous pdf approximaton
    % g should not sum to 1, if the distribution is outside the domain
    
    g = foperator(f, mu(i));
    g = g(1:length(f)) / delta;
    
    trapz(samples, g)
    figure(f2h);
    plot(samples, g); hold on; axis square;
end


%% Create figure handle

legends = cell(length(mu)+1,1);
legends{1} = '$f(x)$';

for i = 1:length(mu)
    if (mu(i) < 0.1)
        legends{i+1} = sprintf('$\\mu = %0.2f$', mu(i));
    elseif (mu(i) < 1.0)
        legends{i+1} = sprintf('$\\mu = %0.1f$', mu(i));
    else
        legends{i+1} = sprintf('$\\mu = %0.0f$', mu(i));
    end
end


%% Print to file

figure(f2h);
axis([0 250 10^-4 10^-1]);
xlabel('$x$','interpreter','latex');
ylabel('$g(x)$','interpreter','latex');

l = legend(legends); legend('boxoff');
set(l, 'interpreter','latex');

filename = 'negbinomdist.pdf';
print_cmd = sprintf('../paperfigs/%s', filename); 
print(f2h, print_cmd, '-dpdf', '-painters');
system(sprintf('pdfcrop --margins 1 ../paperfigs/%s ../paperfigs/%s', filename, filename));

