% calculate_metrics.m
function [mse, psnr] = calculate_metrics(original, reconstructed)
    % Calculates MSE and PSNR [cite: 493-497]
    original = double(original);
    reconstructed = double(rescale(reconstructed, 0, 255));
    
    J = numel(original);
    
    % Mean Squared Error
    mse = sum((original(:) - reconstructed(:)).^2) / J;
    
    % Peak Signal-to-Noise Ratio
    max_I = 255; % As per paper [cite: 496]
    psnr = 10 * log10(max_I^2 / mse);
end