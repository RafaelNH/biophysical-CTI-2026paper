close all
clear all
clc

fig = figure('color', [1 1 1], 'Units', 'centimeters', ...
    'Position', [1 1 21 15]);

set(fig, 'PaperPositionMode', 'auto');
set(fig, 'PaperOrientation', 'portrait');
set(groot, 'DefaultTextFontName', 'Arial');

seed_value = 123;
rng(seed_value, 'twister');

D0 = 2.5;
r0 = 1.0;
As = 0.6;
SD_r = As/sqrt(2);        % Standard deviation of the radius
lp = 1.5 + 4*r0;          % Length of the bead (µm)
rmin = r0 - sqrt(2)*SD_r; % Radius of the connecting tube

% colorcode
num_lines = 6;
c_start = [45, 11, 117]/255;  
c_end = [233, 184, 98]/255; 
colors = [linspace(c_start(1), c_end(1), num_lines)', ...
          linspace(c_start(2), c_end(2), num_lines)', ...
          linspace(c_start(3), c_end(3), num_lines)'];

% Titles for the subplots
titles_up = {{'A1) Periodic', '(Sinusoidal Beads)'},...
             {'A2) Periodic', '(Beads + Tubes)'}, ...
             {'A3) Non-Periodic', '(Disordered)'}};

titles_down = {'B1)','B2)','B3)'};

for cases=1:3
    
    subplot(2, 3, cases)
    
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
    
    z_plot = linspace(0, L_z, 5000); % High res for long super-cell
    idx_plot = discretize(z_plot, z_edges);
    is_bead_plot = mod(idx_plot, 2) == 1;
    
    r_plot = repmat(rmin, size(z_plot));
    z_local_plot = z_plot - z_edges(idx_plot)';
    r_plot(is_bead_plot) = r0 - sqrt(2)*SD_r * cos(2*pi * z_local_plot(is_bead_plot) / lp);
    
    % --- Plot: The Geometry Profile ---
    plot(z_plot, r_plot, 'k-', 'linewidth', 1.5); hold on;
    plot(z_plot, -r_plot, 'k-', 'linewidth', 1.5)
    fill([z_plot, fliplr(z_plot)], [r_plot, fliplr(-r_plot)], [0.9 0.9 0.9], 'EdgeColor', 'none')
    
    %title(sprintf('Disordered Beaded Axon (Variance = %.1f)', var_tube_length), 'Interpreter', 'latex', 'FontSize', 14)
    %xlabel('Axial Position $z$ ($\mu m$)', 'Interpreter', 'latex', 'FontSize', 12)
    %ylabel('Radius ($\mu m$)', 'Interpreter', 'latex', 'FontSize', 12)
    ylim([-max_r*1.2, max_r*1.2])
    xlim([0, min(L_z, lp*10)]) % Zoom in on the first few periods for visibility
    title(titles_up{cases}, 'FontSize', 12)
    grid off; axis equal; hold off;
    
    
    % --- Plot: Axial Kurtosis ---- %
    subplot(2, 3, 3 + cases)
    hold all
    
    legend_str = cell(1, 6);
    color_idx = 1;
    
    for As=0.1:0.1:0.6
        SD_r = As/sqrt(2);
        if cases == 1
        legend_str{color_idx} = sprintf('$$SD(r) = %.2f \\ \\mu m$$', SD_r);
        end
        load(['sims_',num2str(cases),'_As_', num2str(As*100)])
        plot(time_steps_saved/ctime, Kzz_vec, 'color', colors(color_idx, :), 'linewidth', 2)
        
        color_idx = color_idx + 1;
    end
    % Add Y-label to the first column
    if cases == 1
        ylabel('$$K^{\parallel}_{\mu}$$', 'Interpreter', 'latex', 'FontSize', 14)
    end
    xlabel('Diffusion Time  $\Delta / \Delta_c$', 'Interpreter', 'latex', 'FontSize', 12)
    
    % Add the legend
    if cases==1
    leg = legend(legend_str, 'Interpreter', 'latex', 'Location', 'northeast');
    leg.FontSize = 9;
    leg.ItemTokenSize = [15, 18]; % Keeps the colored lines in the legend tidy
    legend('boxoff')
    end
    title(titles_down{cases}, 'FontSize', 12)
    grid on; axis square;
    ylim([0 3.5])
end


print(fig, 'figure7.pdf', '-dpdf', '-painters');
%print(fig, 'figure7.eps', '-depsc', '-painters');
fprintf('Figure saved successfully!\n');