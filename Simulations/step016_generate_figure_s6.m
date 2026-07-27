close all
clear all
clc

load Noise_simulations.mat

fig_width_cm = 25;
fig_height_cm = 20;
fig = figure('color', [1 1 1], 'Units', 'centimeters', ...
    'Position', [-30 -5 fig_width_cm fig_height_cm]);
set(fig, 'PaperPositionMode', 'auto');
set(fig, 'PaperOrientation', 'portrait');
set(groot, 'DefaultTextFontName', 'Arial');
lw = 2;
fff = 10;
fff0 = 12;
fff2 = 10;

%% Mean Kurtosis
haxes1 = subplot(3, 3, 1);
hold on

% 25%, 50% (median), and 75% percentiles across the noise repetitions
MK_p25 = prctile(MK_noise, 25, 2);
MK_median = prctile(MK_noise, 50, 2);
MK_p75 = prctile(MK_noise, 75, 2);

% shaded interquartile range (IQR)
x_fill = [AS(:)'/sqrt(2), fliplr(AS(:)')/sqrt(2)];
y_fill = [MK_p25(:)', fliplr(MK_p75(:)')];

% plot 
fill(haxes1, x_fill, y_fill, [0.8 0.8 0.8], 'EdgeColor', 'none', 'FaceAlpha', 0.5, 'DisplayName', 'IQR (25% - 75%)');
plot(haxes1, AS/sqrt(2), MK_median(:), 'color', [60 167 60]/255, 'LineWidth', lw, 'DisplayName', 'Median');
plot(haxes1, AS/sqrt(2), MK(:), 'k--', 'LineWidth', lw, 'DisplayName', 'Noise-free Reference');
xlim([0, 0.6/sqrt(2)])
ylim([0.8 1.6])
title({'A) Mean Total Kurtosis'}, 'fontsize', fff2)
xlabel('$SD(r)/r_0$', 'Interpreter','latex', 'fontsize', fff)
ylabel('$$\overline{K}_{t}$$', 'Interpreter','latex', 'fontsize', fff)
set(haxes1, 'FontSize', fff);
box(haxes1, 'on');
legend(haxes1, 'Location', 'best', 'FontSize', fff-2);


%% Radial Variance Kurtosis
haxes1 = subplot(3, 3, 2);
hold on

% 25%, 50% (median), and 75% percentiles across the noise repetitions
RK_p25 = prctile(RK_noise, 25, 2);
RK_median = prctile(RK_noise, 50, 2);
RK_p75 = prctile(RK_noise, 75, 2);

% shaded interquartile range (IQR)
x_fill = [AS(:)'/sqrt(2), fliplr(AS(:)')/sqrt(2)];
y_fill = [RK_p25(:)', fliplr(RK_p75(:)')];

% plot 
fill(haxes1, x_fill, y_fill, [0.8 0.8 0.8], 'EdgeColor', 'none', 'FaceAlpha', 0.5, 'DisplayName', 'IQR (25% - 75%)');
plot(haxes1, AS/sqrt(2), RK_median(:), 'color', [30 117 30]/255, 'LineWidth', lw, 'DisplayName', 'Median');
plot(haxes1, AS/sqrt(2), RK(:), 'k--', 'LineWidth', lw, 'DisplayName', 'Noise-free Reference');
xlim([0, 0.6/sqrt(2)])
ylim([2 3.2])
title({'B) Radial Total Kurtosis'}, 'fontsize', fff2)
xlabel('$SD(r)/r_0$', 'Interpreter','latex', 'fontsize', fff)
ylabel('$$K_{t}^\bot$$', 'Interpreter','latex', 'fontsize', fff)
set(haxes1, 'FontSize', fff);
box(haxes1, 'on');
legend(haxes1, 'Location', 'best', 'FontSize', fff-2);


%% Axial Total Kurtosis
haxes1 = subplot(3, 3, 3);
hold on

% 25%, 50% (median), and 75% percentiles across the noise repetitions
AK_p25 = prctile(AK_noise, 25, 2);
AK_median = prctile(AK_noise, 50, 2);
AK_p75 = prctile(AK_noise, 75, 2);

% shaded interquartile range (IQR)
x_fill = [AS(:)'/sqrt(2), fliplr(AS(:)')/sqrt(2)];
y_fill = [AK_p25(:)', fliplr(AK_p75(:)')];

% plot 
fill(haxes1, x_fill, y_fill, [0.8 0.8 0.8], 'EdgeColor', 'none', 'FaceAlpha', 0.5, 'DisplayName', 'IQR (25% - 75%)');
plot(haxes1, AS/sqrt(2), AK_median(:), 'color', [100 217 100]/255, 'LineWidth', lw, 'DisplayName', 'Median');
plot(haxes1, AS/sqrt(2), AK(:), 'k--', 'LineWidth', lw, 'DisplayName', 'Noise-free Reference');
xlim([0, 0.6/sqrt(2)])
ylim([-0.2 0.8])
title({'C) Axial Total Kurtosis'}, 'fontsize', fff2)
xlabel('$SD(r)/r_0$', 'Interpreter','latex', 'fontsize', fff)
ylabel('$$K_{t}^\parallel$$', 'Interpreter','latex', 'fontsize', fff)
set(haxes1, 'FontSize', fff);
box(haxes1, 'on');
legend(haxes1, 'Location', 'best', 'FontSize', fff-2);



%% Mean Variance Kurtosis
haxes1 = subplot(3, 3, 4);
hold on

% 25%, 50% (median), and 75% percentiles across the noise repetitions
MKv_p25 = prctile(MKv_noise, 25, 2);
MKv_median = prctile(MKv_noise, 50, 2);
MKv_p75 = prctile(MKv_noise, 75, 2);

% shaded interquartile range (IQR)
x_fill = [AS(:)'/sqrt(2), fliplr(AS(:)')/sqrt(2)];
y_fill = [MKv_p25(:)', fliplr(MKv_p75(:)')];

% plot 
fill(haxes1, x_fill, y_fill, [0.8 0.8 0.8], 'EdgeColor', 'none', 'FaceAlpha', 0.5, 'DisplayName', 'IQR (25% - 75%)');
plot(haxes1, AS/sqrt(2), MKv_median(:), 'color', [167 82 56]/255, 'LineWidth', lw, 'DisplayName', 'Median');
plot(haxes1, AS/sqrt(2), MKv(:), 'k--', 'LineWidth', lw, 'DisplayName', 'Noise-free Reference');
xlim([0, 0.6/sqrt(2)])
ylim([0.8 1.6])
title({'D) Mean Variance Kurtosis'}, 'fontsize', fff2)
xlabel('$SD(r)/r_0$', 'Interpreter','latex', 'fontsize', fff)
ylabel('$$\overline{K}_{v}$$', 'Interpreter','latex', 'fontsize', fff)
set(haxes1, 'FontSize', fff);
box(haxes1, 'on');
legend(haxes1, 'Location', 'best', 'FontSize', fff-2);


%% Radial Variance Kurtosis
haxes1 = subplot(3, 3, 5);
hold on

% 25%, 50% (median), and 75% percentiles across the noise repetitions
RKv_p25 = prctile(RKv_noise, 25, 2);
RKv_median = prctile(RKv_noise, 50, 2);
RKv_p75 = prctile(RKv_noise, 75, 2);

% shaded interquartile range (IQR)
x_fill = [AS(:)'/sqrt(2), fliplr(AS(:)')/sqrt(2)];
y_fill = [RKv_p25(:)', fliplr(RKv_p75(:)')];

% plot 
fill(haxes1, x_fill, y_fill, [0.8 0.8 0.8], 'EdgeColor', 'none', 'FaceAlpha', 0.5, 'DisplayName', 'IQR (25% - 75%)');
plot(haxes1, AS/sqrt(2), RKv_median(:), 'color', [117 11 45]/255, 'LineWidth', lw, 'DisplayName', 'Median');
plot(haxes1, AS/sqrt(2), RKv(:), 'k--', 'LineWidth', lw, 'DisplayName', 'Noise-free Reference');
xlim([0, 0.6/sqrt(2)])
ylim([2 3.2])
title({'E) Radial Variance Kurtosis'}, 'fontsize', fff2)
xlabel('$SD(r)/r_0$', 'Interpreter','latex', 'fontsize', fff)
ylabel('$$K_{v}^\bot$$', 'Interpreter','latex', 'fontsize', fff)
set(haxes1, 'FontSize', fff);
box(haxes1, 'on');
legend(haxes1, 'Location', 'best', 'FontSize', fff-2);


%% Axial Variance Kurtosis
haxes1 = subplot(3, 3, 6);
hold on

% 25%, 50% (median), and 75% percentiles across the noise repetitions
AKv_p25 = prctile(AKv_noise, 25, 2);
AKv_median = prctile(AKv_noise, 50, 2);
AKv_p75 = prctile(AKv_noise, 75, 2);

% shaded interquartile range (IQR)
x_fill = [AS(:)'/sqrt(2), fliplr(AS(:)')/sqrt(2)];
y_fill = [AKv_p25(:)', fliplr(AKv_p75(:)')];

% plot 
fill(haxes1, x_fill, y_fill, [0.8 0.8 0.8], 'EdgeColor', 'none', 'FaceAlpha', 0.5, 'DisplayName', 'IQR (25% - 75%)');
plot(haxes1, AS/sqrt(2), AKv_median(:), 'color', [217 154 68]/255, 'LineWidth', lw, 'DisplayName', 'Median');
plot(haxes1, AS/sqrt(2), AKv(:), 'k--', 'LineWidth', lw, 'DisplayName', 'Noise-free Reference');
xlim([0, 0.6/sqrt(2)])
ylim([-.2 .8])
title({'F) Axial Variance Kurtosis'}, 'fontsize', fff2)
xlabel('$SD(r)/r_0$', 'Interpreter','latex', 'fontsize', fff)
ylabel('$$K_{v}^\parallel$$', 'Interpreter','latex', 'fontsize', fff)
set(haxes1, 'FontSize', fff);
box(haxes1, 'on');
legend(haxes1, 'Location', 'best', 'FontSize', fff-2);


%% Mean  Microscopic Kurtosis
haxes1 = subplot(3, 3, 7);
hold on

% 25%, 50% (median), and 75% percentiles across the noise repetitions
MKi_p25 = prctile(MKi_noise, 25, 2);
MKi_median = prctile(MKi_noise, 50, 2);
MKi_p75 = prctile(MKi_noise, 75, 2);

% shaded interquartile range (IQR)
x_fill = [AS(:)'/sqrt(2), fliplr(AS(:)')/sqrt(2)];
y_fill = [MKi_p25(:)', fliplr(MKi_p75(:)')];

% plot 
fill(haxes1, x_fill, y_fill, [0.8 0.8 0.8], 'EdgeColor', 'none', 'FaceAlpha', 0.5, 'DisplayName', 'IQR (25% - 75%)');
plot(haxes1, AS/sqrt(2), MKi_median(:), 'color', [56 82 167]/255, 'LineWidth', lw, 'DisplayName', 'Median');
plot(haxes1, AS/sqrt(2), MKi(:), 'k--', 'LineWidth', lw, 'DisplayName', 'Noise-free Reference');
xlim([0, 0.6/sqrt(2)])
ylim([-0.2 0.8])
title({'G) Mean Microscopic Kurtosis'}, 'fontsize', fff2)
xlabel('$SD(r)/r_0$', 'Interpreter','latex', 'fontsize', fff)
ylabel('$$\overline{K}_{\mu}$$', 'Interpreter','latex', 'fontsize', fff)
set(haxes1, 'FontSize', fff);
box(haxes1, 'on');
legend(haxes1, 'Location', 'best', 'FontSize', fff-2);


%% Radial  Microscopic Kurtosis
haxes1 = subplot(3, 3, 8);
hold on

% 25%, 50% (median), and 75% percentiles across the noise repetitions
RKi_p25 = prctile(RKi_noise, 25, 2);
RKi_median = prctile(RKi_noise, 50, 2);
RKi_p75 = prctile(RKi_noise, 75, 2);

% shaded interquartile range (IQR)
x_fill = [AS(:)'/sqrt(2), fliplr(AS(:)')/sqrt(2)];
y_fill = [RKi_p25(:)', fliplr(RKi_p75(:)')];

% plot 
fill(haxes1, x_fill, y_fill, [0.8 0.8 0.8], 'EdgeColor', 'none', 'FaceAlpha', 0.5, 'DisplayName', 'IQR (25% - 75%)');
plot(haxes1, AS/sqrt(2), RKi_median(:), 'color', [45 11 117]/255, 'LineWidth', lw, 'DisplayName', 'Median');
plot(haxes1, AS/sqrt(2), RKi(:), 'k--', 'LineWidth', lw, 'DisplayName', 'Noise-free Reference');
xlim([0, 0.6/sqrt(2)])
ylim([-0.2 0.8])
title({'H) Radial Microscopic Kurtosis'}, 'fontsize', fff2)
xlabel('$SD(r)/r_0$', 'Interpreter','latex', 'fontsize', fff)
ylabel('$$K_{\mu}^\bot$$', 'Interpreter','latex', 'fontsize', fff)
set(haxes1, 'FontSize', fff);
box(haxes1, 'on');
legend(haxes1, 'Location', 'north', 'FontSize', fff-2);


%% Axial Microscopic Kurtosis
haxes1 = subplot(3, 3, 9);
hold on

% 25%, 50% (median), and 75% percentiles across the noise repetitions
AKi_p25 = prctile(AKi_noise, 25, 2);
AKi_median = prctile(AKi_noise, 50, 2);
AKi_p75 = prctile(AKi_noise, 75, 2);

% shaded interquartile range (IQR)
x_fill = [AS(:)'/sqrt(2), fliplr(AS(:)')/sqrt(2)];
y_fill = [AKi_p25(:)', fliplr(AKi_p75(:)')];

% plot 
fill(haxes1, x_fill, y_fill, [0.8 0.8 0.8], 'EdgeColor', 'none', 'FaceAlpha', 0.5, 'DisplayName', 'IQR (25% - 75%)');
plot(haxes1, AS/sqrt(2), AKi_median(:), 'color', [68 154 217]/255, 'LineWidth', lw, 'DisplayName', 'Median');
plot(haxes1, AS/sqrt(2), AKi(:), 'k--', 'LineWidth', lw, 'DisplayName', 'Noise-free Reference');
xlim([0, 0.6/sqrt(2)])
ylim([-0.2 0.8])
title({'I) Axial Microscopic Kurtosis'}, 'fontsize', fff2)
xlabel('$SD(r)/r_0$', 'Interpreter','latex', 'fontsize', fff)
ylabel('$$K_{\mu}^\parallel$$', 'Interpreter','latex', 'fontsize', fff)
set(haxes1, 'FontSize', fff);
box(haxes1, 'on');
legend(haxes1, 'Location', 'best', 'FontSize', fff-2);

print(fig, 'figureS6.pdf', '-dpdf', '-painters');
fprintf('Figures saved successfully!\n');