close all
clear all
clc

% Note to run these code you have to generate the inside and outside
% kurtosis contributions from step005

fs = filesep;
addpath(['..', fs, 'DWIu_toolbox_v1.11', fs, 'DKI_toolbox/'])
addpath(['..', fs, 'DWIu_toolbox_v1.11', fs, 'DDE_toolbox/'])
addpath(['..', fs, 'DWIu_toolbox_v1.11', fs, 'Spherical_functions_and_directions/'])
load('dirs1024.mat')

lw = 1.5;
ls = '--';
fff = 10;
fff2 = 10;
fig = figure('color', [1 1 1], 'Units', 'centimeters', ...
    'Position', [1 1 21 8]);

set(fig, 'PaperPositionMode', 'auto');
set(fig, 'PaperOrientation', 'portrait');
set(groot, 'DefaultTextFontName', 'Arial');

%% DDE profiles
www = 0.24;
for si = 1:2
    handle = subplot(1, 2, si);
    %handle = subplot('position', [0.04+(si-1)*www, 0.15, 0.20, 0.7]);
    if si == 1
        load('sim_kurt_vs_as_inside', 'AS',  'RKi', 'AKi', 'RKv', 'AKv', 'RK', 'AK')
        plot(AS/sqrt(2), RKi,'LineWidth', lw, 'color', [45 11 117]/255)
        hold on
        plot(AS/sqrt(2), AKi,'LineWidth', lw, 'color', [68 154 217]/255)
        yline(0, 'LineWidth', lw, 'color', [0.8 0.8 0.8])
        xticks(0:0.1:0.4)
        yticks(-0.1:0.2:0.5)
        xlim([0, AS(end)/sqrt(2)])
        ylim([-0.1, 0.5])        
        legend('$$K_{\mu}^\bot$$', '$$K_{\mu}^\parallel$$', '$$y=0$$', ...
            'Interpreter','latex', 'fontsize', fff, 'location', 'northwest')
        title('A) Inside Contribution', 'fontsize', fff2)
        xlabel('$SV(r) (\mu m)$', 'Interpreter','latex', 'fontsize', fff)
        
    elseif si==2
        load('sim_kurt_vs_as_outside', 'AS',  'RKi', 'AKi', 'RKv', 'AKv', 'RK', 'AK')
        plot(AS/sqrt(2), RKi,'LineWidth', lw, 'color', [45 11 117]/255)
        hold on
        plot(AS/sqrt(2), AKi,'LineWidth', lw, 'color', [68 154 217]/255)
        yline(0, 'LineWidth', lw, 'color', [0.8 0.8 0.8])
        xticks(0:0.1:0.4)
        yticks(-0.1:0.2:0.5)
        xlim([0, AS(end)/sqrt(2)])
        ylim([-0.1, 0.5])       
        legend('$$K_{\mu}^\bot$$', '$$K_{\mu}^\parallel$$', '$$y=0$$', ...
            'Interpreter','latex', 'fontsize', fff, 'location', 'northwest')
        title('B) Outside Contribution', 'fontsize', fff2)
        xlabel('$SV(r) (\mu m)$', 'Interpreter','latex', 'fontsize', fff)
    end
end

print(fig, 'figureS4.pdf', '-dpdf', '-painters');
%print(fig, 'figureS1_final.eps', '-depsc', '-painters');
fprintf('Figures saved successfully!\n');
