% main_experiment.m
% Main script to reimplement the numerical experiment from:
% "Regularization in a functional reproducing kernel Hilbert space"
% by R. Wang and Y. Xu, Journal of Complexity 66 (2021) 101567.

clear; clc; close all;
rng('default'); % for reproducibility

%% --- 1. Experiment Setup ---
disp('Setting up experiment...');

% Image and Radon Transform Parameters
img_size = 512; % Paper uses 512, using 256 for speed.
m = 4096;       % Number of Radon projections. Paper uses 16,384.
phantom_name = 'shepp-logan'; % 'bullseye', 'crescent', 'bubbles', 'shepp-logan'

% Generate grid for image reconstruction
[X, Y] = meshgrid(linspace(-1, 1, img_size), linspace(-1, 1, img_size));
pixel_grid = [X(:), Y(:)];

% Generate Radon transform coordinates (randomly, as per the paper)
% MATLAB's radon function needs angles in degrees
theta_deg = 180 * rand(m, 1);
z_k = [cosd(theta_deg), sind(theta_deg)]; % Projection directions
t_k = (2 * rand(m, 1) - 1); % Projection positions in [-1, 1] relative to image diagonal

%% --- 2. Data Generation ---
disp('Generating data...');

% Get the original phantom image
original_image = generate_phantoms(phantom_name, img_size);
figure; imshow(original_image, []); title('Original Phantom');

% Generate noise-free Radon data using MATLAB's radon function
disp('Generating sinogram...');
% The 't_k' positions need to be scaled to the pixel dimensions for the radon function
[R, xp] = radon(original_image, theta_deg);
y_noise_free = zeros(m, 1);
for i = 1:m
    % Interpolate to get the radon value at our specific t_k
    y_noise_free(i) = interp1(xp, R(:,i), t_k(i) * max(xp), 'linear', 0);
end

% Generate noisy data
noise_std = 0.02 * (max(y_noise_free) - min(y_noise_free));
noise = noise_std * randn(m, 1);
y_noisy = y_noise_free + noise;

%% --- 3. Reconstruction ---
disp('Starting reconstruction...');

% Method parameters are taken from Tables 7 & 8 for the Shepp-Logan phantom
% --- FR Method ---
disp('Running FR method...');
params_FR.lambda = 3.2258e-5; % From Table 7 (noise-free)
params_FR.delta = 3.5484;
params_FR.gamma = 0.0007;
img_FR = reconstruct_FR(y_noise_free, z_k, t_k, pixel_grid, img_size, params_FR);
[mse_FR, psnr_FR] = calculate_metrics(original_image, img_FR);
fprintf('FR (Noise-Free): MSE = %.4e, PSNR = %.2f dB\n', mse_FR, psnr_FR);

% --- RR Method ---
disp('Running RR method...');
params_RR.lambda = 1.0e-6; % From Table 7
params_RR.gamma = 0.0020;
img_RR = reconstruct_RR(y_noise_free, z_k, t_k, pixel_grid, img_size, params_RR);
[mse_RR, psnr_RR] = calculate_metrics(original_image, img_RR);
fprintf('RR (Noise-Free): MSE = %.4e, PSNR = %.2f dB\n', mse_RR, psnr_RR);

% --- LR-G Method ---
disp('Running LR-G method...');
params_LRG.lambda = 9.6865; % From Table 7
params_LRG.gamma = 0.0010;
img_LRG = reconstruct_LR_G(y_noise_free, z_k, t_k, pixel_grid, img_size, params_LRG);
[mse_LRG, psnr_LRG] = calculate_metrics(original_image, img_LRG);
fprintf('LR-G (Noise-Free): MSE = %.4e, PSNR = %.2f dB\n', mse_LRG, psnr_LRG);


%% --- 4. Display Results ---
figure;
subplot(2, 2, 1); imshow(original_image, []); title('Original');
subplot(2, 2, 2); imshow(img_FR, []); title(sprintf('FR (PSNR: %.2f dB)', psnr_FR));
subplot(2, 2, 3); imshow(img_RR, []); title(sprintf('RR (PSNR: %.2f dB)', psnr_RR));
subplot(2, 2, 4); imshow(img_LRG, []); title(sprintf('LR-G (PSNR: %.2f dB)', psnr_LRG));
sgtitle('Reconstructions from Noise-Free Data');
