%close all; clear; clc;

fs = filesep;
addpath(['..', fs, 'DWIu_toolbox_v1.11', fs, 'Spherical_functions_and_directions/'])

% formating figure page and fonts
fig2_width_cm = 25;
fig2_height_cm = 25;
%fig2 = figure('color', [1 1 1]);
fig2 = figure('color', [1 1 1], 'Units', 'centimeters', ...
    'Position', [-30 -5 fig2_width_cm fig2_height_cm]);
set(fig2, 'PaperPositionMode', 'auto');
set(fig2, 'PaperOrientation', 'portrait');
set(groot, 'DefaultTextFontName', 'Arial');
lw = 1.5;
fff = 10;
fff0 = 12;


% Plot ODF for different OP values
haxes0 = subplot(3, 1, 1);
haxes0.Position = [0.1300 0.65, 1*0.7, 0.3*0.7];
np = 200;
theta = linspace(0,2*pi,np+1);
phi = linspace(0,2*pi,np+1);
[theta, phi] = meshgrid(theta, phi);

% ODF plot
OPs = [0.1:0.2:0.9];
kappas = OP2kappa(OPs);

k_ind = 0;
for k=kappas;
k_ind = k_ind + 1;
angle = asin(sqrt(1/(2*k)))/pi*180;

B = Bingham_2D_for_ploting(1, k, k, [0, 1, 0]', [1, 0, 0]', theta, phi);
[xx, yy, zz] = sph2cart(theta, phi, B);

surf(xx+2*(k_ind-3), yy, zz, B, 'EdgeColor', 'none')
hold on
view(0,0)

text(2*(k_ind-3)-1, 1, 1, ['OP=',num2str(OPs(k_ind))],...
    'Interpreter','latex', 'fontsize', fff0)
end

title({'A) ODFs for different OP values'}, 'fontsize', 12)
zlim([-1.5, 1.5])
xlim([-5, 5])
ylim([-5, 5])
axis off
camlight   
lighting gouraud



% Results of sensitivity to OP in the beaded case
load('sim_kurt_vs_op_bead')

haxes1 = subplot(3, 4, 5);
plot(OPs, RKv,'LineWidth', lw, 'color', [117 11 45]/255)
hold on
plot(OPs, AKv,'LineWidth', lw, 'color', [217 154 68]/255)
plot(OPs, kaniso+kiso,'LineWidth', lw, 'color', [0 0 0],'Linestyle', ':')
xticks([0:0.2:1])
yticks([0:1:4])
xlim([0.01, 0.9])
ylim([0, 5])
legend('$$K_{v}^\bot$$', '$$K_{v}^\parallel$$', ...
            'Interpreter','latex', 'fontsize', fff, 'location', 'north')
xlabel('$OP$', 'Interpreter','latex', 'fontsize', fff)
text(0.05, 4.5, 'B1', 'fontsize', fff)
text(0.4, kaniso(1)+kiso(1)+0.3, '$$\overline{K}_{ani}+\overline{K}_{iso}$$',...
    'Interpreter','latex', 'fontsize', fff0)

subplot(3, 4, 6)
plot(OPs, RKi,'LineWidth', lw, 'color', [45 11 117]/255)
hold on
plot(OPs, AKi,'LineWidth', lw, 'color', [68 154 217]/255)
plot(OPs, MKi,'LineWidth', lw, 'color', [0 0 0],'Linestyle', ':')
xticks([0:0.2:1])
yticks([0:0.2:1])
xlim([0.01, 0.9])
ylim([0, 0.7])
        legend('$$K_{\mu}^\bot$$', '$$K_{\mu}^\parallel$$', ...
            'Interpreter','latex', 'fontsize', fff, 'location', 'north')
xlabel('$OP$', 'Interpreter','latex', 'fontsize', fff)
text(0.05, 4.5/5*0.7, 'B2', 'fontsize', fff)
text(0.75, kintra(1)-0.4/5*0.7, '$$\overline{K}_{\mu}$$',...
    'Interpreter','latex', 'fontsize', fff0)
title({'B) Sensitivity to dispersion (SD(r)/r_0=0.42, No Exchange)'}, 'fontsize', fff)%, 'HorizontalAlignment', 'center')

subplot(3, 4, 7)
RK = RKv + RKi;
AK = AKv + AKi;
plot(OPs, RKv./RK,'LineWidth', lw, 'color', [117 11 45]/255)
hold on
plot(OPs, AKv./AK,'LineWidth', lw, 'color', [217 154 68]/255, 'Linestyle', ':')
plot(OPs, RKi./RK,'LineWidth', lw, 'color', [45 11 117]/255)
plot(OPs, AKi./AK,'LineWidth', lw, 'color', [68 154 217]/255, 'Linestyle', ':')
xticks([0:0.2:1])
yticks([0:0.2:1])
xlim([0.01, 0.9])
ylim([0, 1.1])
%legend('$$K_{v}^\bot/K_t^\bot$$', '$$K_{v}^\parallel/K_t^\parallel$$', ...
%            '$$K_{\mu}^\bot/K_t^\bot$$', '$$K_{\mu}^\parallel/K_t^\parallel$$', ...
%            'Interpreter','latex', 'fontsize', fff, 'location', 'east')
xlabel('$OP$', 'Interpreter','latex', 'fontsize', fff)
text(0.05, 4.5/5*1.1, 'B3', 'fontsize', fff)

subplot(3, 4, 8)
% Dummy subplot for legend
plot(OPs, RKv./RK,'LineWidth', lw, 'color', [117 11 45]/255)
hold on
plot(OPs, AKv./AK,'LineWidth', lw, 'color', [217 154 68]/255, 'Linestyle', ':')
plot(OPs, RKi./RK,'LineWidth', lw, 'color', [45 11 117]/255)
plot(OPs, AKi./AK,'LineWidth', lw, 'color', [68 154 217]/255, 'Linestyle', ':')
%xticks([0:0.2:1])
%yticks([0:0.2:1])
xlim([0.01, 0.9])
ylim([10, 20])
axis off
lgd = legend('$$K_{v}^\bot/K_t^\bot$$', '$$K_{v}^\parallel/K_t^\parallel$$', ...
            '$$K_{\mu}^\bot/K_t^\bot$$', '$$K_{\mu}^\parallel/K_t^\parallel$$', ...
            'Interpreter','latex', 'fontsize', fff, 'location', 'west');
lgd.Position(1) = lgd.Position(1) - 0.045;
%xlabel('$OP$', 'Interpreter','latex', 'fontsize', fff)
%text(0.05, 4.5/5*1.1, 'B3', 'fontsize', fff)


% Results of sensitivity to OP in the exchange case
load('sim_kurt_vs_op_exchange')

subplot(3, 4, 9)
plot(OPs, RKv,'LineWidth', lw, 'color', [117 11 45]/255)
hold on
plot(OPs, AKv,'LineWidth', lw, 'color', [217 154 68]/255)
plot(OPs, kaniso+kiso,'LineWidth', lw, 'color', [0 0 0],'Linestyle', ':')
xticks([0:0.2:1])
yticks([0:1:4])
xlim([0.01, 0.9])
ylim([0, 2.2])
legend('$$K_{v}^\bot$$', '$$K_{v}^\parallel$$', ...
            'Interpreter','latex', 'fontsize', fff, 'location', 'north')
xlabel('$OP$', 'Interpreter','latex', 'fontsize', fff)
text(0.05, 4.5/5*2.2, 'C1', 'fontsize', fff)
text(0.4, kaniso(1)+kiso(1)-0.04, '$$\overline{K}_{ani}+\overline{K}_{iso}$$',...
    'Interpreter','latex', 'fontsize', fff0)

subplot(3, 4, 10)
plot(OPs, RKi,'LineWidth', lw, 'color', [45 11 117]/255)
hold on
plot(OPs, AKi,'LineWidth', lw, 'color', [68 154 217]/255)
plot(OPs, kintra,'LineWidth', lw, 'color', [0 0 0],'Linestyle', ':')
xticks([0:0.2:1])
yticks([0:1:4])
xlim([0.01, 0.9])
ylim([0, 2.2])
        legend('$$K_{\mu}^\bot$$', '$$K_{\mu}^\parallel$$', ...
            'Interpreter','latex', 'fontsize', fff, 'location', 'north')
xlabel('$OP$', 'Interpreter','latex', 'fontsize', fff)
text(0.05, 4.5/5*2.2, 'C2', 'fontsize', fff)
text(0.75, kintra(end)+0.35/5*2.2, '$$\overline{K}_{\mu}$$',...
    'Interpreter','latex', 'fontsize', fff0)
title({'C) Sensitivity to dispersion (SD(r)/r_0=0, k=0.3/\Delta)'}, 'fontsize', fff)%, 'HorizontalAlignment', 'center')

subplot(3, 4, 11)
RK = RKv + RKi;
AK = AKv + AKi;
plot(OPs, RKv./RK,'LineWidth', lw, 'color', [117 11 45]/255)
hold on
plot(OPs, AKv./AK,'LineWidth', lw, 'color', [217 154 68]/255, 'Linestyle', ':')
plot(OPs, RKi./RK,'LineWidth', lw, 'color', [45 11 117]/255)
plot(OPs, AKi./AK,'LineWidth', lw, 'color', [68 154 217]/255, 'Linestyle', ':')
xticks([0:0.2:1])
yticks([0:0.2:1])
xlim([0.01, 0.9])
ylim([0, 1.1])
%legend('$$K_{v}^\bot/K_t^\bot$$', '$$K_{v}^\parallel/K_t^\parallel$$', ...
%            '$$K_{\mu}^\bot/K_t^\bot$$', '$$K_{\mu}^\parallel/K_t^\parallel$$', ...
%            'Interpreter','latex', 'fontsize', fff, 'location', 'east')
xlabel('$OP$', 'Interpreter','latex', 'fontsize', fff)
text(0.05, 4.5/5*1.1, 'C3', 'fontsize', fff)

subplot(3, 4, 12)
% Dummy subplot for legend
plot(OPs, RKv./RK,'LineWidth', lw, 'color', [117 11 45]/255)
hold on
plot(OPs, AKv./AK,'LineWidth', lw, 'color', [217 154 68]/255, 'Linestyle', ':')
plot(OPs, RKi./RK,'LineWidth', lw, 'color', [45 11 117]/255)
plot(OPs, AKi./AK,'LineWidth', lw, 'color', [68 154 217]/255, 'Linestyle', ':')
%xticks([0:0.2:1])
%yticks([0:0.2:1])
xlim([0.01, 0.9])
ylim([10, 20])
axis off
lgd = legend('$$K_{v}^\bot/K_t^\bot$$', '$$K_{v}^\parallel/K_t^\parallel$$', ...
            '$$K_{\mu}^\bot/K_t^\bot$$', '$$K_{\mu}^\parallel/K_t^\parallel$$', ...
            'Interpreter','latex', 'fontsize', fff, 'location', 'west');
lgd.Position(1) = lgd.Position(1) - 0.045;
%xlabel('$OP$', 'Interpreter','latex', 'fontsize', fff)
%text(0.05, 4.5/5*1.1, 'B3', 'fontsize', fff)

print(fig2, 'figure2_final.pdf', '-dpdf', '-painters');
%print(fig2, 'figure2_final.eps', '-depsc', '-painters');
fprintf('Figures saved successfully!\n');
