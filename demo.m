%% run_demo.m
% Demo for LtFCTN tensor completion


clc; clear; close all;
addpath(genpath(pwd));
rng('default');
rng(2,'twister');

%% ================= Load data =================
load('color_suzie.mat')
X = X(:,:,:,1:50);
X = double(X);
if max(X(:))>1
    X = X/max(X(:));
end

Nway = size(X);

%% ================= Settings =================
sample_ratio = 0.20;

fprintf('\n=============================================\n');
fprintf(' LtFCTN demo for tensor completion\n');
fprintf(' Tensor size: ');
fprintf('%d ', Nway);
fprintf('\n Sampling ratio: %.2f\n', sample_ratio);
fprintf('=============================================\n\n');

%% ================= Generate mask =================
Omega = find(rand(prod(Nway),1) < sample_ratio);
Y = zeros(Nway);
Y(Omega) = X(Omega);
Mask = false(Nway);
Mask(Omega) = true;
Mask_perm = permute(Mask, [1,2,4,3]);
Omega_LtFCTN = find(Mask_perm);

%% ================= Observed quality =================
[psnr_obs, ssim_obs, fsim_obs] = my_quality(X*255, Y*255);

%% ================= LtFCTN parameters =================
opts = [];
opts.P = [75, 80, 22, 2];

Tp  = opts.P(end);
r12 = 28;
r13 = 5;
r14 = 7;

opts.R = [0,   r12, r13, Tp;
          0,   0,   r14, Tp;
          0,   0,   0,   Tp;
          Tp,  Tp,  Tp,  Tp];

opts.tol   = 1e-5;
opts.maxit = 500;
opts.rho   = 0.001;
opts.mu    = 10;
opts.Xtrue = X;

%% ================= Run LtFCTN =================
fprintf('Running LtFCTN...\n');

t0 = tic;
[X_rec, G, Out] = LtFCTN(Y, Omega_LtFCTN, opts);
time_LtFCTN = toc(t0);

X_rec = min(max(X_rec,0),1);

%% ================= Evaluation =================
[psnr_LtFCTN, ssim_LtFCTN, fsim_LtFCTN] = my_quality(X*255, X_rec*255);

fprintf('\n================== Results ==================\n');
fprintf('%12s | %8s | %8s | %8s | %8s\n', ...
    'Method','PSNR','SSIM','FSIM','Time(s)');
fprintf('%12s | %8.4f | %8.4f | %8.4f | %8s\n', ...
    'Observed', psnr_obs, ssim_obs, fsim_obs, '-');
fprintf('%12s | %8.4f | %8.4f | %8.4f | %8.2f\n', ...
    'LtFCTN', psnr_LtFCTN, ssim_LtFCTN, fsim_LtFCTN, time_LtFCTN);
fprintf('=============================================\n');


%% ================= Save results =================
if ~exist('results','dir')
    mkdir('results');
end

save(fullfile('results','LtFCTN_demo_result.mat'), ...
    'X_rec', 'G', 'Out', ...
    'sample_ratio', ...
    'psnr_obs', 'ssim_obs', 'fsim_obs', ...
    'psnr_LtFCTN', 'ssim_LtFCTN', 'fsim_LtFCTN', ...
    'time_LtFCTN');

fprintf('\nResults saved to results/LtFCTN_demo_result.mat\n');