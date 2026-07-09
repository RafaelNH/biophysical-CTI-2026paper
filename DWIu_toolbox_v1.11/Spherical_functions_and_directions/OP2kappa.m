function kappa = OP2kappa(OP)

    % Initialize kappa array
    kappa = zeros(size(OP));
    OP_input = OP;
    OP_input(OP_input < 0) = 0;
    OP_input(OP_input > 1) = 1; 

    % Forward Lookup Table (LUT)
    k_grid = linspace(0, 700, 10000)';
    x = linspace(0, 1, 2000); 
    integrand = exp(k_grid * (x.^2)); 
    I0_grid = trapz(x, integrand, 2);
    
    cos2_theta_grid = (exp(k_grid) ./ (2 .* k_grid .* I0_grid)) - (1 ./ (2 .* k_grid));
    OP_grid = (3 .* cos2_theta_grid - 1) ./ 2;
    OP_grid(1) = 0; % Fix mathematical singularity for k=0
    
    % Ensure the grid is monotonic for inverse interpolation
    [OP_uniq, uniq_idx] = unique(OP_grid);
    k_uniq = k_grid(uniq_idx);

    % Interpolation (0 <= OP < 0.99)
    valid_idx = (OP_input >= 0) & (OP_input < 0.99);
    if any(valid_idx, 'all')
        kappa(valid_idx) = interp1(OP_uniq, k_uniq, OP_input(valid_idx), 'pchip');
    end
    
    % Invert extreme cases
    large_idx = (OP_input >= 0.99) & (OP_input < 1);
    if any(large_idx, 'all')
        kappa(large_idx) = 3 ./ (4 .* (1 - OP_input(large_idx)));
    end

    % infinite concentration for OP == 1
    inf_idx = (OP_input == 1);
    if any(inf_idx, 'all')
        kappa(inf_idx) = inf;
    end
end