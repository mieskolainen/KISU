% Simple Inversion of ALICE data: <https://www.hepdata.net/record/ins1614477>
%
% mikael.mieskolainen@cern.ch, 2019

clear; close all;
addpath('./src_hist');
addpath('./src');

system('mkdir ../lhcfigs');
system('mkdir ../lhcfigs2D');

%% Parameter setup

% Poisson mean hypothesis values
mu_vals = [0.33 1.5 3.5 6];

param.algo   = 5;    % Algorithm 5 for Poisson case, 3 for Gaussian approx
param.R      = 100;  % Number of iterations in the inverse algorithm (>= 100)
param.lambda = 0.03; % Fixed regularization strength

N_bootstrap  = 100;   % Number of bootstrap samples

%{
%param.principle = 'a_minus_ndf_min'; % Discrepancy principle
param.principle = 'a_plus_b_min';     % Variational equilibrium
%param.principle = 'a_minus_b_min';   % Test

% Lambda values to scan
lambdas  = logspace(-2.301, 0.1, 30);
%}

% Skip 0-bin
SKIP0BIN = false;

% Domain extension (technical factor, at least 2 x to be safe)
extension = 2.0;

%% Read data

data_read;

%% Data inverse

data_invert;

%% Data inverse 2D

data_invert_2D;

%% mu-value scans

data_mu_scan;
