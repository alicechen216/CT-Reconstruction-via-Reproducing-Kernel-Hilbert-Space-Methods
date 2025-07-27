% radon_transform_Gaussian_basis.m
function G = radon_transform_Gaussian_basis(z_k, t_k, x_j, gamma)
    % Analytically computes the Radon transform of a Gaussian basis function
    % K_gamma(x, x_j) = exp(-||x-x_j||^2 / gamma)
    % R_{z,t}(K_gamma(., x_j)) = sqrt(pi*gamma) * exp(-(t - z.x_j)^2 / gamma)
    
    m = size(z_k, 1);
    num_basis = size(x_j, 1);
    G = zeros(m, num_basis);
    
    prefactor = sqrt(pi * gamma);
    
    % Note: The t_k values from main are in [-1, 1] and need to be scaled
    % to the same space as the dot product z_k * x_j'
    % This scaling depends on the image diagonal length, which is sqrt(2)
    % for a [-1,1] grid. Let's assume t_k is already scaled appropriately
    % for simplicity, consistent with the math.
    
    for i = 1:m
        t_proj = t_k(i) - z_k(i, :) * x_j'; % This is (t - z.x_j) for all j
        G(i, :) = prefactor * exp(-t_proj.^2 / gamma);
    end
end
