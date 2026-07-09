function OP = kappa2OP(kappa)

    % Initialize OP array
    OP = zeros(size(kappa));

    % For valid numerical computation (0 < kappa < 700)
    valid_idx = (kappa > 0) & (kappa < 700);
    k_input = kappa(valid_idx);
    
    if ~isempty(k_input)
        % Lookup Table
        k_grid = linspace(0, 700, 10000)';
        
        % integral for normalization
        x = linspace(0, 1, 2000);
        integrand = exp(k_grid * (x.^2)); 
        I0_grid = trapz(x, integrand, 2);
        
        % Calculate cos^2(theta) for the grid
        cos2_theta_grid = (exp(k_grid) ./ (2 .* k_grid .* I0_grid)) - (1 ./ (2 .* k_grid));
        
        % Calculate OP for the grid
        OP_grid = (3 .* cos2_theta_grid - 1) ./ 2;
        
        % Fix the mathematical singularity exactly at k = 0
        OP_grid(1) = 0; 

        % Interpolate OP against the LUT
        OP(valid_idx) = interp1(k_grid, OP_grid, k_input, 'linear', 'extrap');
    end
    
    %  Cases that kappa will introduce overflow
    large_idx = (kappa >= 700);
    if any(large_idx, 'all')
        OP(large_idx) = 1;
    end
end