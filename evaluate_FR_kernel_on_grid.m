% evaluate_FR_kernel_on_grid.m
function K_R_grid = evaluate_FR_kernel_on_grid(z_k, t_k, pixel_grid, params)
    % Evaluates K_R(z_k, t_k)(x) for all k and all pixels x in pixel_grid.
    % K_R from Eq. (4.24)
    delta = params.delta;
    gamma = params.gamma;
    m = size(z_k, 1);
    num_pixels = size(pixel_grid, 1);
    K_R_grid = zeros(num_pixels, m);

    A = delta + 1/gamma;
    C = 1 / (gamma * (1 + gamma * delta));

    prefactor = sqrt(pi / A);

    for i = 1:m
        zk = z_k(i, :);
        tk = t_k(i);
        zk_perp = [-zk(2), zk(1)];

        x_dot_z = pixel_grid * zk';
        x_dot_z_perp = pixel_grid * zk_perp';
        norm_x_sq = sum(pixel_grid.^2, 2);

        exponent = -A * (norm_x_sq + tk^2) + (2*tk/gamma) * x_dot_z + C * x_dot_z_perp.^2;
        K_R_grid(:, i) = prefactor * exp(exponent);
    end
end
