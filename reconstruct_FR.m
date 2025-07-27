% reconstruct_FR.m
function g_lambda = reconstruct_FR(y, z_k, t_k, pixel_grid, img_size, params)
    % Reconstructs an image using the FR method described in the paper.
    % Solves (F_R + lambda*I)c = y
    % g_lambda = sum(c_j * K_R(z_j, t_j))

    lambda = params.lambda;
    m = length(y);

    % --- Step 1: Assemble the matrix F_R ---
    % F_R(j, k) = R_{z_j, t_j}(K_R(z_k, t_k, .))
    fprintf('Assembling FR matrix...\n');
    F_R = radon_transform_FR_kernel(z_k, t_k, z_k, t_k, params);

    % --- Step 2: Solve the linear system for coefficients c ---
    fprintf('Solving FR linear system...\n');
    I = speye(m);
    A = F_R + lambda * I;
    [c, ~] = pcg(A, y, 1e-4, 200); % Using Conjugate Gradient as per paper

    % --- Step 3: Reconstruct the image on the pixel grid ---
    fprintf('Reconstructing final FR image...\n');
    % g_lambda(x) = sum_j c_j * K_R((z_j, t_j), x)
    K_R_image = evaluate_FR_kernel_on_grid(z_k, t_k, pixel_grid, params);
    g_lambda_vec = K_R_image * c;
    g_lambda = flipud(reshape(g_lambda_vec, img_size, img_size));
end

