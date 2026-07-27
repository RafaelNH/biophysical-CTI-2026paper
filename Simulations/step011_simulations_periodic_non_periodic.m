close all
clear all
clc

% Set the pseudo-random seed for reproducibility
seed_value = 123;
rng(seed_value, 'twister');

disp('Initializing 3D Non-Periodic (Disordered) Beaded Cylinder Simulation...');

%% 1. Simulation Parameters
D0 = 2.5;         % Intrinsic Diffusivity (µm^2/ms)
dt = 1e-3;        % Time step (ms)
npar = 100000;    % Number of particles

% Geometry Parameters
r0 = 1.0;         % Mean radius (µm)
for As = 0.2:-0.1:0.1
    SD_r = As/sqrt(2);      % Standard deviation of the radius
    lp = 1.5 + 4*r0;  % Length of the bead (µm)
    rmin = r0 - sqrt(2)*SD_r; % Radius of the connecting tube
    
    % cases 1 (no tubes),
    % cases 2 (uniform tubes),
    % cases 3 (short range disorder),
    cases = 3;
    if cases == 1
        mean_tube_length = 0;
        var_tube_length  = 0;
        ctime = lp^2/(2*D0); % homogenization time
    elseif cases == 2
        mean_tube_length = lp;
        var_tube_length  = 0;
        ctime = (2*lp)^2/(2*D0);
    else
        mean_tube_length = lp;
        var_tube_length  = 4.0;
        ctime = (2*lp)^2/(2*D0);
    end
    T = ctime * 3;
    
    num_units = 50;   % Size of the "Super-Cell" (number of beads/tubes)
    
    % Generate tube lengths
    if var_tube_length == 0
        tube_lengths = ones(num_units, 1) * mean_tube_length;
    else
        % Gamma distribution parameters: shape (k) and scale (theta)
        shape_k = (mean_tube_length^2) / var_tube_length;
        scale_theta = var_tube_length / mean_tube_length;
        tube_lengths = gamrnd(shape_k, scale_theta, [num_units, 1]);
    end
    
    % Build the 1D map of the Super-Cell edges
    % z_edges will look like: [0, bead1_end, tube1_end, bead2_end, tube2_end...]
    z_edges = zeros(num_units * 2 + 1, 1);
    current_z = 0;
    for i = 1:num_units
        current_z = current_z + lp;
        z_edges(2*i) = current_z; % End of bead
        
        current_z = current_z + tube_lengths(i);
        z_edges(2*i + 1) = current_z; % End of tube
    end
    
    L_z = z_edges(end); % Total length of the Super-Cell
    max_r = r0 + sqrt(2)*SD_r;
    
    % Pre-calculate step standard deviation
    step_std = sqrt(2 * D0 * dt);
    display = false;
    
    %% 2. Initialize Particles
    pos = zeros(npar, 3);
    valid = false(npar, 1);
    disp('Placing particles...');
    
    % --- Plot 1: The Geometry Profile ---
    figure('color', [1 1 1], 'Position', [50, 100, 1200, 350]);
    z_plot = linspace(0, L_z, 5000); % High res for long super-cell
    idx_plot = discretize(z_plot, z_edges);
    is_bead_plot = mod(idx_plot, 2) == 1;
    
    r_plot = repmat(rmin, size(z_plot));
    z_local_plot = z_plot - z_edges(idx_plot)';
    r_plot(is_bead_plot) = r0 - sqrt(2)*SD_r * cos(2*pi * z_local_plot(is_bead_plot) / lp);
    
    plot(z_plot, r_plot, 'k-', 'linewidth', 1.5); hold on;
    plot(z_plot, -r_plot, 'k-', 'linewidth', 1.5)
    fill([z_plot, fliplr(z_plot)], [r_plot, fliplr(-r_plot)], [0.9 0.9 0.9], 'EdgeColor', 'none')
    
    title(sprintf('Disordered Beaded Axon (Variance = %.1f)', var_tube_length), 'Interpreter', 'latex', 'FontSize', 14)
    xlabel('Axial Position $z$ ($\mu m$)', 'Interpreter', 'latex', 'FontSize', 12)
    ylabel('Radius ($\mu m$)', 'Interpreter', 'latex', 'FontSize', 12)
    ylim([-max_r*1.2, max_r*1.2])
    xlim([0, min(L_z, lp*15)]) % Zoom in on the first few periods for visibility
    grid off; axis equal; hold off;
    
    while ~all(valid)
        num_invalid = sum(~valid);
        
        x_rand = (rand(num_invalid, 1) * 2 - 1) * max_r;
        y_rand = (rand(num_invalid, 1) * 2 - 1) * max_r;
        z_rand = rand(num_invalid, 1) * L_z;
        
        idx_rand = discretize(z_rand, z_edges);
        is_bead_rand = mod(idx_rand, 2) == 1;
        z_local_rand = z_rand - z_edges(idx_rand);
        
        rz_local = repmat(rmin, num_invalid, 1);
        rz_local(is_bead_rand) = r0 - sqrt(2)*SD_r * cos(2*pi * z_local_rand(is_bead_rand) / lp);
        
        is_inside = (x_rand.^2 + y_rand.^2) <= rz_local.^2;
        
        temp_pos = [x_rand, y_rand, z_rand];
        pos(~valid, :) = temp_pos;
        valid(~valid) = is_inside;
    end
    
    pos_abs = pos;
    pos0_abs = pos;
    
    %% 3. Pre-allocate arrays
    time_steps = dt:dt:T;
    num_steps = length(time_steps);
    
    % Define sampling frequency
    sample_freq = 10;
    
    % Create an array of only the steps we actually save
    saved_step_indices = sample_freq:sample_freq:num_steps;
    time_steps_saved = time_steps(saved_step_indices);
    num_saved = length(time_steps_saved);
    
    Dzz_vec = zeros(1, num_saved);
    Kzz_vec = zeros(1, num_saved);
    Dxx_vec = zeros(1, num_saved);
    Kxx_vec = zeros(1, num_saved);
    
    save_idx = 1; % Counter to track our position in the saved arrays
    
    disp('Running Random Walk...');
    
    %% 4. Main Simulation Loop
    for ai = 1:num_steps
        st = randn(npar, 3) * step_std;
        pos_trial = pos + st;
        
        % Fast locate particles in the Super-Cell
        z_mod = mod(pos_trial(:,3), L_z);
        idx = discretize(z_mod, z_edges);
        is_bead = mod(idx, 2) == 1;
        z_local = z_mod - z_edges(idx);
        
        % Compute local radii
        rz_trial = repmat(rmin, npar, 1);
        rz_trial(is_bead) = r0 - sqrt(2)*SD_r * cos(2*pi * z_local(is_bead) / lp);
        
        % Collision detection
        r2_trial = pos_trial(:,1).^2 + pos_trial(:,2).^2;
        out_idx = r2_trial > rz_trial.^2;
        
        if any(out_idx)
            z_out_local = z_local(out_idx);
            rz_out = rz_trial(out_idx);
            is_bead_out = is_bead(out_idx);
            
            dr_dz = zeros(sum(out_idx), 1);
            dr_dz(is_bead_out) = sqrt(2)*SD_r * (2*pi/lp) .* sin(2*pi * z_out_local(is_bead_out) / lp);
            
            nx = 2 * pos_trial(out_idx, 1);
            ny = 2 * pos_trial(out_idx, 2);
            nz = -2 * rz_out .* dr_dz;
            
            nn = sqrt(nx.^2 + ny.^2 + nz.^2);
            nx = nx ./ nn;  ny = ny ./ nn;  nz = nz ./ nn;
            
            st_out = st(out_idx, :);
            dot_product = st_out(:,1).*nx + st_out(:,2).*ny + st_out(:,3).*nz;
            
            st_ref_x = st_out(:,1) - 2 * dot_product .* nx;
            st_ref_y = st_out(:,2) - 2 * dot_product .* ny;
            st_ref_z = st_out(:,3) - 2 * dot_product .* nz;
            
            pos_trial(out_idx, 1) = pos(out_idx, 1) + st_ref_x;
            pos_trial(out_idx, 2) = pos(out_idx, 2) + st_ref_y;
            pos_trial(out_idx, 3) = pos(out_idx, 3) + st_ref_z;
            
            % Fallback safety
            z_mod_check = mod(pos_trial(out_idx, 3), L_z);
            idx_check = discretize(z_mod_check, z_edges);
            is_bead_check = mod(idx_check, 2) == 1;
            z_local_check = z_mod_check - z_edges(idx_check);
            
            rz_check = repmat(rmin, sum(out_idx), 1);
            rz_check(is_bead_check) = r0 - sqrt(2)*SD_r * cos(2*pi * z_local_check(is_bead_check) / lp);
            
            still_out = (pos_trial(out_idx,1).^2 + pos_trial(out_idx,2).^2) > rz_check.^2;
            
            if any(still_out)
                actual_out_indices = find(out_idx);
                fail_idx = actual_out_indices(still_out);
                pos_trial(fail_idx, :) = pos(fail_idx, :);
            end
        end
        
        % --- ABSOLUTE TRACKING & PBC WRAP ---
        actual_st = pos_trial - pos;
        pos_abs = pos_abs + actual_st;
        pos_trial(:, 3) = mod(pos_trial(:, 3), L_z);
        pos = pos_trial;
        
        % --- METRICS ---
        if mod(ai, sample_freq) == 0
            current_time = time_steps(ai);
            
            dZ = pos_abs(:, 3) - pos0_abs(:, 3);
            mZ2 = mean(dZ.^2);
            mZ4 = mean(dZ.^4);
            
            Dzz_vec(save_idx) = mZ2 / (2 * current_time);
            Kzz_vec(save_idx) = (mZ4 / (mZ2^2)) - 3;
            
            dX = pos_abs(:, 1) - pos0_abs(:, 1);
            mX2 = mean(dX.^2);
            mX4 = mean(dX.^4);
            Dxx_tmp = mX2 / (2 * current_time);
            Kxx_tmp = (mX4 / (mX2^2)) - 3;
            
            dY = pos_abs(:, 2) - pos0_abs(:, 2);
            mY2 = mean(dY.^2);
            mY4 = mean(dY.^4);
            Dyy_tmp = mY2 / (2 * current_time);
            Kyy_tmp = (mY4 / (mY2^2)) - 3;
            
            Dxx_vec(save_idx) = (Dxx_tmp + Dyy_tmp) / 2;
            Kxx_vec(save_idx) = (Kxx_tmp + Kyy_tmp) / 2;
            
            save_idx = save_idx + 1; % Advance the save index
        end
        
        % (Keep your display and progress tracker code here)
        if display && mod(ai, 50) == 0
            plot(pos(:,3), pos(:,1), '.b', 'MarkerSize', 1)
            ylim([-max_r*1.5, max_r*1.5])
            xlim([0, min(L_z, lp*15)])
            axis equal
            drawnow
        end
        
        if mod(ai, 1000) == 0
            fprintf('Progress: %.1f ms / %.1f ms\n', current_time, T);
        end
    end
    
    disp('Simulation Complete!');
    
    %% 5. Plotting Results
    figure('color', [1 1 1], 'Position', [100, 100, 600, 600]);
    
    subplot(2, 2, 1)
    plot(time_steps_saved / ctime, Dzz_vec, 'color', [68 154 217]/255, 'linewidth', 2)
    title('Apparent Axial Diffusivity ($D_{zz}$)', 'Interpreter', 'latex', 'FontSize', 14)
    xlabel('Diffusion Time  $\Delta$ /$\Delta_c$', 'Interpreter', 'latex', 'FontSize', 12)
    ylabel('$D_{zz}$ ($\mu m^2$/ms)', 'Interpreter', 'latex', 'FontSize', 12)
    grid on; axis square;
    
    subplot(2, 2, 2)
    plot(time_steps_saved/ctime, Kzz_vec, 'color', [45 11 117]/255, 'linewidth', 2)
    title('Axial Micro-Kurtosis ($K_{zz}$)', 'Interpreter', 'latex', 'FontSize', 14)
    xlabel('Diffusion Time  $\Delta$ /$\Delta_c$', 'Interpreter', 'latex', 'FontSize', 12)
    ylabel('$K_{zz}$', 'Interpreter', 'latex', 'FontSize', 12)
    grid on; axis square;
    
    subplot(2, 2, 3)
    plot(time_steps_saved/ctime, Dxx_vec, 'color', [68 154 217]/255, 'linewidth', 2)
    title('Apparent Radial Diffusivity ($D_{xx}$)', 'Interpreter', 'latex', 'FontSize', 14)
    xlabel('Diffusion Time $\Delta$ /$\Delta_c$', 'Interpreter', 'latex', 'FontSize', 12)
    ylabel('$D_{zz}$ ($\mu m^2$/ms)', 'Interpreter', 'latex', 'FontSize', 12)
    grid on; axis square;
    
    subplot(2, 2, 4)
    plot(time_steps_saved/ctime, Kxx_vec, 'color', [45 11 117]/255, 'linewidth', 2)
    title('Radial Micro-Kurtosis ($K_{xx}$)', 'Interpreter', 'latex', 'FontSize', 14)
    xlabel('Diffusion Time  $\Delta$ /$\Delta_c$', 'Interpreter', 'latex', 'FontSize', 12)
    ylabel('$K_{zz}$', 'Interpreter', 'latex', 'FontSize', 12)
    grid on; axis square;
    
    save(['sims_',num2str(cases),'_As_', num2str(As*100)], 'Dzz_vec', 'Kzz_vec', 'Dxx_vec', 'Kxx_vec', 'time_steps_saved', 'ctime')
    disp(['Save sims_',num2str(cases),'_As_', num2str(As*100)])
    close all
    
end