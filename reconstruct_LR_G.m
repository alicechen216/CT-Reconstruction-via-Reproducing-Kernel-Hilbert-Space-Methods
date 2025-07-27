% reconstruct_LR_G.m
function g_lambda = reconstruct_LR_G(y, z_k, t_k, pixel_grid, img_size, params)
    % Reconstructs an image using the LR-G method.
    % Solves (G*G' + lambda*E)c = G*y
    
    lambda = params.lambda;
    gamma = params.gamma;

    % Define the basis centers (uniform grid)
    grid_res = 32; % Paper uses 64, using 32 for speed
    [X_basis, Y_basis] = meshgrid(linspace(-1, 1, grid_res), linspace(-1, 1, grid_res));
    x_j = [X_basis(:), Y_basis(:)];

    % --- Step 1: Assemble matrix G ---
    % Same as in RR method
    fprintf('Assembling LR-G matrix G...\n');
    G = radon_transform_Gaussian_basis(z_k, t_k, x_j, gamma);

    % --- Step 2: Assemble matrix E ---
    % E(j,k) = <phi_k, phi_j>_L2
    % Integral of K_gamma(x,x_j)*K_gamma(x,x_k) over R^2
    % This integral is (pi*gamma/2) * exp(-||x_j-x_k||^2 / (2*gamma))
    fprintf('Assembling LR-G matrix E...\n');
    dist_sq = pdist2(x_j, x_j, 'squaredeuclidean');
    E = (pi * gamma / 2) * exp(-dist_sq / (2 * gamma));

    % --- Step 3: Solve the linear system ---
    fprintf('Solving LR-G linear system...\n');
    A = G' * G + lambda * E;
    b = G' * y;
    [c, ~] = pcg(A, b, 1e-4, 200);

    % --- Step 4: Reconstruct the image ---
    fprintf('Reconstructing final LR-G image...\n');
    K_gamma_image = exp(-pdist2(pixel_grid, x_j, 'squaredeuclidean') / gamma);
    g_lambda_vec = K_gamma_image * c;
    g_lambda = flipud(reshape(g_lambda_vec, img_size, img_size));
end