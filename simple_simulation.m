% Fixed K-fold inverse autoconvolution
% via naive analytic spectral domain inverse versus recursive KISU solution
%
% mikael.mieskolainen@cern.ch, 2019

clear; close all;
addpath('./src_hist');
addpath('./src');

system('mkdir ../paperfigs');

%% Simulation setup

K = 2; % How many times the convolution (CONSTANT == no stochastic)

% True distribution is negative exponential
b = 2;
x_vals = linspace(0.1,20,100);
f = 1/b*exp(-x_vals/b);
f = f/sum(f);

% Do the autoconvolution directly
g = f;
for i = 1:K-1
    g = conv(g,f);
end
g(g < 0 ) = 0;

% Implement the right total number of events
N = 1e3;
g = g / sum(g);
g = g * N;

% Add noise to the measurement (just Gaussian for simplicity)
sigma = 0.5;
g = g + randn(size(g))*sigma;

% Make sure it is non-negative
g(g < 0 ) = 0;


%% Naive FFT solution with high frequency hard cutoff

CUTOFF = [0.4];
fhat = {};
for i = 1:length(CUTOFF)
    fhat{i} = naivefft(g, x_vals, K, CUTOFF(i));
end


%% KISU iterative solution

mu = zeros(K,1); mu(end) = 1; % Construct fixed mu vector

lambda = 0.1;  % Regularization
ITER   = 100;  % Iterations
algo   = 5;    % Algorithm

[f_iters, rho, F_tot, g_err, reg_err, g_hat] = kisusolver(g, [], mu, ITER, lambda, algo);
f_hat_iter = f_iters(1:length(f),end);


%% Plot

fig1 = figure;

% x-values
xx = linspace(x_vals(1),x_vals(end)*K, length(g));

for pl = 2:-1:1
    subplot(1,2,pl);
    plot(x_vals, f* N, 'k-.'); hold on;
    stepbar(xx, g, 'r-');
    
    for i = 1:length(fhat)
        stepbar(x_vals, fhat{i}, '-', 'color', [0.6 0.6 0.6]*i, 'linewidth', 1);
    end
    
    stepbar(x_vals, f_hat_iter, '-', 'color', [0 0 0.9]);
    
    if (pl == 1)
        axis([0 x_vals(end)*0.6 0 N*0.1]);
        ylabel('counts','interpreter','latex');
    end
    if (pl == 2) 
        set(gca,'yscale','log');
        axis([0 x_vals(end)*0.6 1 N*0.1]);
    end
    axis square;
    xlabel('$x$','interpreter','latex');
    set(gca,'XTick', 0:2:x_vals(end));
end

legs = {'True $f$', sprintf('Smeared $g$, $K = %d$', K)};
for i = 1:length(fhat)
   legs{end+1} = sprintf('Naive $\\Lambda_\\omega = %0.2f$', CUTOFF(i));
end
legs{end+1} = sprintf('KISU $\\lambda = %0.2f$', lambda);

l = legend(legs);
set(l, 'interpreter','latex'); legend('boxoff');

filename = sprintf('simplesimu1.pdf');
print_cmd = sprintf('print -dpdf ../paperfigs/%s', filename); 
eval(print_cmd);
system(sprintf('pdfcrop --margins 1 ../paperfigs/%s ../paperfigs/%s', filename, filename));


%% Look at the Fourier transform

F = fft(f);

figure;
subplot(2,2,1);
stem(real(F));
xlabel('Re[F]','interpreter','latex'); axis square;

subplot(2,2,2);
stem(imag(F));
xlabel('Im[F]','interpreter','latex'); axis square;

subplot(2,2,3);
stem(abs(F));
xlabel('$|F|$','interpreter','latex'); axis square;

subplot(2,2,4);
stem(phase(F));
xlabel('arg[F]','interpreter','latex'); axis square;
%}
