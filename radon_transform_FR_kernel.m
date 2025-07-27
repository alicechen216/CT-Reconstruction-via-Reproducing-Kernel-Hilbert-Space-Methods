% radon_transform_FR_kernel.m
function F_R = radon_transform_FR_kernel(z_rows, t_rows, z_cols, t_cols, params)
    % Analytically computes R_{z_j, t_j}(K_R(z_k, t_k, .))
    % This is the matrix F_R where F_R(j,k) is the transform.
    % The derivation involves solving a Gaussian integral.
    
    delta = params.delta;
    gamma = params.gamma;
    
    num_rows = size(z_rows, 1);
    num_cols = size(z_cols, 1);
    F_R = zeros(num_rows, num_cols);

    A = delta + 1/gamma;
    C = 1 / (gamma * (1 + gamma * delta));
    
    prefactor_k = sqrt(pi / A);

    for j = 1:num_rows % Index for the transform R_{z_j, t_j}
        zj = z_rows(j, :);
        tj = t_rows(j);
        zj_perp = [-zj(2), zj(1)];
        
        for k = 1:num_cols % Index for the kernel K_R(z_k, t_k)
            zk = z_cols(k, :);
            tk = t_cols(k);
            zk_perp = [-zk(2), zk(1)];
            
            % Terms for the quadratic in 's' (the integration variable)
            % Exponent is P*s^2 + Q*s + R
            P = -A + C * dot(zj_perp, zk_perp)^2 + (2/gamma) * dot(zj_perp, zk) * dot(zj_perp, zk_perp) * C;
            P = -A + C * dot(zj_perp, zk_perp)^2; % Simplified from paper structure

            Q = 2*tj*C*dot(zj, zk_perp)*dot(zj_perp, zk_perp) + ...
                (2/gamma)*(tk*dot(zj_perp, zk) + tj*dot(zj,zk)*dot(zj_perp,zk));

            Q = 2*tj*C*dot(zj, zk_perp)*dot(zj_perp, zk_perp) + (2*tk/gamma)*dot(zj_perp, zk);

            R = -A*(tj^2 + tk^2) + (2*tk/gamma)*tj*dot(zj, zk) + C*(tj*dot(zj, zk_perp))^2;

            % The integral of exp(P*s^2 + Q*s + R) is sqrt(-pi/P) * exp(R - Q^2/(4P))
            integral_val = sqrt(pi / -P) * exp(R - (Q^2 / (4 * P)));
            
            F_R(j, k) = prefactor_k * integral_val;
        end
    end
end
