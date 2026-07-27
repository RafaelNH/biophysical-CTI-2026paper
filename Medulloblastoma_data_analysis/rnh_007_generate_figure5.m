clear all
close all
clc

fig = figure('color', [1 1 1], 'Units', 'centimeters', ...
    'Position', [-30 -5 21 23]);
set(fig, 'PaperPositionMode', 'auto');
set(fig, 'PaperOrientation', 'portrait');
set(groot, 'DefaultTextFontName', 'Arial');

fff = 12;


%% Histology
annotation('textbox', [0.05, 0.81, 0.4, 0.1], ...
            'String', '$$\textbf{A) Histology (Control)}$$', ...
            'HorizontalAlignment', 'center', ...
            'BackgroundColor', 'none', ...
            'EdgeColor', 'none', 'Interpreter', 'latex', 'fontsize', fff);

annotation('textbox', [0.52, 0.81, 0.4, 0.1], ...
            'String', '$$\textbf{B) Histology (Medulloblastoma)}$$', ...
            'HorizontalAlignment', 'center', ...
            'BackgroundColor', 'none', ...
            'EdgeColor', 'none', 'Interpreter', 'latex', 'fontsize', fff);

fs = filesep;

P1 = [0.04, 0.65, 0.20, 0.28];
P2 = [0.28, 0.65, 0.20, 0.28];
P3 = [0.52, 0.65, 0.20, 0.28];
P4 = [0.76, 0.65, 0.20, 0.28];

ax1 = subplot_tight(3, 4, 1);

Image2 = imread('panelA1.jpg');
imagesc(Image2)
hold on
wi = 175;
xi = 290;
yi = 82;
plot([xi, xi, xi+wi, xi+wi, xi],[yi yi+wi yi+wi yi yi], 'color', 'black')
axis equal
axis off
text(120, 90, 'A1)', 'BackGroundColor', [0 0 0], 'Interpreter','latex', 'fontsize', 10, 'color', [1 1 1], 'HorizontalAlignment', 'right')
text(800, 810, '$$2mm$$', 'BackGroundColor', 'none', 'Interpreter','latex', 'fontsize', 10, 'color', [0 0 0], 'HorizontalAlignment', 'right')
ax1.Position = P1;

ax1 = subplot_tight(3, 4, 2);
Image2 = imread('panelA2.jpg');
imagesc(Image2)
axis equal
axis off
text(140, 90, 'A2)', 'BackGroundColor', [0 0 0], 'Interpreter','latex', 'fontsize', 10, 'color', [1 1 1], 'HorizontalAlignment', 'right')
text(800, 810, '$$500{\mu}m$$', 'BackGroundColor', 'none', 'Interpreter','latex', 'fontsize', 10, 'color', [0 0 0], 'HorizontalAlignment', 'right')
ax1.Position = P2;

ax1 = subplot_tight(3, 4, 3);
Image2 = imread('panelB1.jpg');
imagesc(Image2)
hold on
xi = 380;
yi = 55;
plot([xi, xi, xi+wi, xi+wi, xi],[yi yi+wi yi+wi yi yi], 'color', 'black')
axis equal
axis off
text(120, 90, 'B1)', 'BackGroundColor', [0 0 0], 'Interpreter','latex', 'fontsize', 10, 'color', [1 1 1], 'HorizontalAlignment', 'right')
text(800, 810, '$$2mm$$', 'BackGroundColor', 'none', 'Interpreter','latex', 'fontsize', 10, 'color', [0 0 0], 'HorizontalAlignment', 'right')
ax1.Position = P3;

ax1 = subplot_tight(3, 4, 4);
Image2 = imread('panelB2.jpg');
imagesc(Image2)
axis equal
axis off
text(140, 90, 'B2)', 'BackGroundColor', [0 0 0], 'Interpreter','latex', 'fontsize', 10, 'color', [1 1 1], 'HorizontalAlignment', 'right')
text(800, 810, '$$500{\mu}m$$', 'BackGroundColor', 'none', 'Interpreter','latex', 'fontsize', 10, 'color', [0 0 0], 'HorizontalAlignment', 'right')
ax1.Position = P4;

%% Analysis 
rnh_006_paper_analysis
fff = 12;
fff2 = 10;

%% Panels C
mycmap = turbo(1000);

%FA
ax = subplot_tight(6, 6, 13);
imagesc(squeeze(FAcontrol), [0 1]); axis image; axis off;
hold on
%title('C1) $$FA$$', 'Interpreter','latex', 'fontsize', fff)
colormap_set
cb.Ticks = [0 1];
text(40, 85, 'C1)  $$FA$$', 'BackGroundColor', 'none', 'Interpreter','latex', 'fontsize', fff, 'color', [1 1 1], 'HorizontalAlignment', 'center')


% KANISO
ax = subplot_tight(6, 6, 14);
imagesc(squeeze(KANIcontrol), [0 1.5]); axis image; axis off;
hold on
%title({'C2) $$\overline{K}_{ani}$$'}, 'Interpreter','latex', 'fontsize', fff)
%title({'\textbf{C) CTI maps (control)}'}, 'Interpreter','latex', 'fontsize', fff, 'HorizontalAlignment', 'left'),...
colormap_set
cb.Ticks = [0 1.5];
text(40, 85, 'C2) $$\overline{K}_{ani}$$', 'BackGroundColor', 'none', 'Interpreter','latex', 'fontsize', fff, 'color', [1 1 1], 'HorizontalAlignment', 'center')
%title({'\textbf{C) CTI maps (control)}'}, 'Interpreter','latex', 'fontsize', fff, 'HorizontalAlignment', 'l'),...


% Kiso
ax = subplot_tight(6, 6, 15);
imagesc(squeeze(KISOcontrol), [0 1.5]); axis image; axis off;
hold on
%title({'\textbf{C) CTI maps (control)}'; 'C3) $$\overline{K}_{iso}$$'}, 'Interpreter','latex', 'fontsize', fff, 'HorizontalAlignment', 'left')
%title({'\textbf{C) CTI maps (control)}'}, 'Interpreter','latex', 'fontsize', fff, 'HorizontalAlignment', 'left'),...
  % 'VerticalAlignment', 'bottom')
colormap_set
cb.Ticks = [0 1.5];
text(40, 85, 'C3) $$\overline{K}_{iso}$$', 'BackGroundColor', 'none', 'Interpreter','latex', 'fontsize', fff, 'color', [1 1 1], 'HorizontalAlignment', 'center')


% Ku
ax = subplot_tight(6, 6, 16);
imagesc(squeeze(MKicontrol), [0 1.5]); axis image; axis off;
hold on
%title('C4) $$\overline{K}_{\mu}$$', 'Interpreter','latex', 'fontsize', fff)
colormap_set
cb.Ticks = [0 1.5];
text(40, 85, 'C4) $$\overline{K}_{\mu}$$', 'BackGroundColor', 'none', 'Interpreter','latex', 'fontsize', fff, 'color', [1 1 1], 'HorizontalAlignment', 'center')



% RKi
ax = subplot_tight(6, 6, 17);
imagesc(squeeze(RKicontrol), [0 1.5]); axis image; axis off;
hold on
%title('C5) $$K^\bot_{\mu}$$', 'Interpreter','latex', 'fontsize', fff)
colormap_set
cb.Ticks = [0 1.5];
text(40, 85, 'C5) $$K^\bot_{\mu}$$', 'BackGroundColor', 'none', 'Interpreter','latex', 'fontsize', fff, 'color', [1 1 1], 'HorizontalAlignment', 'center')


% AKi
ax = subplot_tight(6, 6, 18);
imagesc(squeeze(AKicontrol), [0 1.5]); axis image; axis off;
hold on
%title('C6) $$K^\parallel_{\mu}$$', 'Interpreter','latex', 'fontsize', fff)
colormap_set
cb.Ticks = [0 1.5];
text(40, 85, 'C6) $$K^\parallel_{\mu}$$', 'BackGroundColor', 'none', 'Interpreter','latex', 'fontsize', fff, 'color', [1 1 1], 'HorizontalAlignment', 'center')

% Overarching Title C
annotation('textbox', [0.05, 0.62, 0.90, 0.05], ...
            'String', '$$\textbf{C) CTI maps (Control)}$$', ...
            'HorizontalAlignment', 'center', ...
            'BackgroundColor', 'none', ...
            'EdgeColor', 'none', 'Interpreter', 'latex', 'fontsize', fff);


%% Panels D
%FA
ax = subplot_tight(6, 6, 19);
imagesc(squeeze(FAtumor), [0 1]); axis image; axis off;
hold on
%title('D1) $$FA$$', 'Interpreter','latex', 'fontsize', fff)
colormap_set
cb.Ticks = [0 1];
text(40, 85, 'D1) $$FA$$', 'BackGroundColor', 'none', 'Interpreter','latex', 'fontsize', fff, 'color', [1 1 1], 'HorizontalAlignment', 'center')

% KANISO
ax = subplot_tight(6, 6, 20);
imagesc(squeeze(KANItumor), [0 1.5]); axis image; axis off;
hold on
%title({'D2) $$\overline{K}_{ani}$$'}, 'Interpreter','latex', 'fontsize', fff)
colormap_set
cb.Ticks = [0 1.5];
text(40, 85, 'D2) $$\overline{K}_{ani}$$', 'BackGroundColor', 'none', 'Interpreter','latex', 'fontsize', fff, 'color', [1 1 1], 'HorizontalAlignment', 'center')

% Kiso
ax = subplot_tight(6, 6, 21);
imagesc(squeeze(KISOtumor), [0 1.5]); axis image; axis off;
hold on
%title({'\textbf{D) CTI maps (Medulloblastoma)}'}, 'Interpreter','latex', 'fontsize', fff, 'HorizontalAlignment', 'l')
colormap_set
cb.Ticks = [0 1.5];
text(40, 85, 'D3) $$\overline{K}_{iso}$$', 'BackGroundColor', 'none', 'Interpreter','latex', 'fontsize', fff, 'color', [1 1 1], 'HorizontalAlignment', 'center')

% Ku
ax = subplot_tight(6, 6, 22);
imagesc(squeeze(MKitumor), [0 1.5]); axis image; axis off;
hold on
%title('D4) $$\overline{K}_{\mu}$$', 'Interpreter','latex', 'fontsize', fff)
%title({'\textbf{D) CTI maps (Medulloblastoma)}'}, 'Interpreter','latex', 'fontsize', fff, 'HorizontalAlignment', 'c')
colormap_set
cb.Ticks = [0 1.5];
text(40, 85, 'D4) $$\overline{K}_{\mu}$$', 'BackGroundColor', 'none', 'Interpreter','latex', 'fontsize', fff, 'color', [1 1 1], 'HorizontalAlignment', 'center')

%title({'\textbf{D) CTI maps (Medulloblastoma)}'}, 'Interpreter','latex', 'fontsize', fff, 'HorizontalAlignment', 'r')


% RKi
ax = subplot_tight(6, 6, 23);
imagesc(squeeze(RKitumor), [0 1.5]); axis image; axis off;
hold on
%title('D5) $$K^\bot_{\mu}$$', 'Interpreter','latex', 'fontsize', fff)
colormap_set
cb.Ticks = [0 1.5];
text(40, 85, 'D5) $$K^\bot_{\mu}$$', 'BackGroundColor', 'none', 'Interpreter','latex', 'fontsize', fff, 'color', [1 1 1], 'HorizontalAlignment', 'center')



% AKi
ax = subplot_tight(6, 6, 24);
imagesc(squeeze(AKitumor), [0 1.5]); axis image; axis off;
hold on
%title('D6) $$K^\parallel_{\mu}$$', 'Interpreter','latex', 'fontsize', fff)
colormap_set
cb.Ticks = [0 1.5];
text(40, 85, 'D6) $$K^\parallel_{\mu}$$', 'BackGroundColor', 'none', 'Interpreter','latex', 'fontsize', fff, 'color', [1 1 1], 'HorizontalAlignment', 'center')


% Overarching Title D
annotation('textbox', [0.05, 0.46, 0.90, 0.05], ...
            'String', '$$\textbf{D) CTI maps (Medulloblastoma)}$$', ...
            'HorizontalAlignment', 'center', ...
            'BackgroundColor', 'none', ...
            'EdgeColor', 'none', 'Interpreter', 'latex', 'fontsize', fff);

%%
P1 = [0.04, 0.04, 0.20, 0.24];
P2 = [0.28, 0.04, 0.20, 0.24];
P3 = [0.52, 0.04, 0.20, 0.24];
P4 = [0.76, 0.04, 0.20, 0.24];

%% bar analysis
set(groot,'defaultAxesTickLabelInterpreter','latex');

c = 8;
x = 1/(2*c+5);
xs = (c+1)*x/2;
fff = 10;

ct_color = [0.5 0.5 0.5];
st_color = [1 1 1];

bp = 1:4;
sub_control = [ones(1, 6); 2*ones(1, 6); 3*ones(1, 6); 4*ones(1, 6)];
sub_tumor = [ones(1, 5); 2*ones(1, 5); 3*ones(1, 5); 4*ones(1, 5)];

%% E1) plot mean kurtosis quatities
hAxes = subplot(3, 4, 9);
hAxes.Position = P1;

controls = [2, 5, 6, 8, 9, 11];
tumor = [1, 3, 4, 7, 10];

ct_all = [mk(controls(:)); kaniso(controls(:)); kiso(controls(:)); uk(controls(:))];
st_all = [mk(tumor(:)); kaniso(tumor(:)); kiso(tumor(:)); uk(tumor(:))];
ct = mean(ct_all, 2);
st = mean(st_all, 2);

hAxes.TickLabelInterpreter = 'latex';
b1 = bar(bp-xs, ct, c*x);
b1.FaceColor = ct_color;
hold on
b2 = bar(bp+xs, st, c*x);
b2.FaceColor = st_color;
plot([bp(:)-xs, bp(:)+xs]', [ct, st]', 'black')
plot(sub_control(:)-xs, ct_all(:), '.black')
plot(sub_tumor(:)+xs, st_all(:), '.black')
xlim([0.5, 4.5])
ylim([0, 2])
set(gca,'XTick', [1 2 3 4], 'YTick',[0 0.5 1 1.5 2 2.5 3]);
xticks([1 2 3 4])
xticklabels({'$$\overline{K}_t$$', ...
    '$$\overline{K}_{ani}$$',...
    '$$\overline{K}_{iso}$$',...
    '$$\overline{K}_{\mu}$$'})
title('E1) Mean Quantities', 'Interpreter', 'latex', ...
    'fontsize', fff2)
legend('Controls', 'Medulloblastomas', 'fontsize', 8)
hAxes.FontSize = fff2;

% Stats
% --- MK Metric ---
[h, ptotal, ci, stats] = ttest2(mk(controls(:)), mk(tumor(:)));
cohens_d_mk = (mean(mk(controls(:))) - mean(mk(tumor(:)))) / stats.sd;
disp(['MK: ', num2str(cohens_d_mk)]);

% --- K_aniso Metric ---
[h, pvar1, ci, stats] = ttest2(kaniso(controls(:)), kaniso(tumor(:)));
cohens_d_kaniso = (mean(kaniso(controls(:))) - mean(kaniso(tumor(:)))) / stats.sd;
disp(['K_aniso: ', num2str(cohens_d_kaniso)]);

% --- K_iso Metric ---
[h, pvar2, ci, stats] = ttest2(kiso(controls(:)), kiso(tumor(:)));
cohens_d_kiso = (mean(kiso(controls(:))) - mean(kiso(tumor(:)))) / stats.sd;
disp(['K_iso: ', num2str(cohens_d_kiso)]);

% --- uK Metric ---
[h, pmicro, ci, stats] = ttest2(uk(controls(:)), uk(tumor(:)));
cohens_d_uk = (mean(uk(controls(:))) - mean(uk(tumor(:)))) / stats.sd;
disp(['uK: ', num2str(cohens_d_uk)]);

% pvalues
disp(ptotal)
disp(pvar1)
disp(pvar2)
disp(pmicro)


adj_p = [ptotal, pvar1, pvar2, pmicro]*4;
raw_pvals = [ptotal, pvar1, pvar2, pmicro];

[h_fdr, crit_p, adj_ci_cvrg, adj_p] = fdr_bh(raw_pvals, 0.05, 'pdep', 'no');

for mi = 1:4
    if mi == 1
        y_max = max(mk(:)) + 0.1; % Adjust 0.1 as needed for spacing
    elseif mi == 2
        y_max = max(kaniso(:)) + 0.1;
    elseif mi == 3
        y_max = max(kiso(:)) + 0.1;
    elseif mi == 4
        y_max = max(uk(:)) + 0.1;
    end
    
    if adj_p(mi) < 0.001
        sig_marker = '***';
    elseif adj_p(mi) < 0.01
        sig_marker = '**';
    elseif adj_p(mi) < 0.05
        sig_marker = '*';
    else
        sig_marker = '';
    end
    
    if ~isempty(sig_marker)
        text(mi, y_max, sig_marker, ...
            'HorizontalAlignment', 'center', ...
            'VerticalAlignment', 'bottom', ...
            'FontSize', fff, ...
            'Interpreter', 'latex', ...
            'FontWeight', 'bold')
    end
end



%% E2) plot radial kurtosis quatities

c = 6;
x = 1/(2*c+5);
xs = (c+1)*x/2;

bp = 1:3;
sub_control = [ones(1, 6); 2*ones(1, 6); 3*ones(1, 6)];
sub_tumor = [ones(1, 5); 2*ones(1, 5); 3*ones(1, 5)];

hAxes = subplot_tight(3, 4, 10);
hAxes.Position = P2;

ct_all = [rk(controls(:)); vrk(controls(:)); urk(controls(:))];
st_all = [rk(tumor(:)); vrk(tumor(:)); urk(tumor(:))];
ct = mean(ct_all, 2);
st = mean(st_all, 2);

hAxes.TickLabelInterpreter = 'latex';
b1 = bar(bp-xs, ct, c*x);
b1.FaceColor = ct_color;
hold on
b2 = bar(bp+xs, st, c*x);
b2.FaceColor = st_color;
plot([bp(:)-xs, bp(:)+xs]', [ct, st]', 'black')
plot(sub_control(:)-xs, ct_all(:), '.black')
plot(sub_tumor(:)+xs, st_all(:), '.black')
xlim([0.5, 3.5])
ylim([0, 2])
set(gca,'XTick', [1 2 3], 'YTick',[0 0.5 1 1.5 2 2.5 3]);
xticks([1 2 3])
xticklabels({'$$K^\bot_t$$', ...
    '$$K^\bot_v$$',...
    '$$K^\bot_{\mu}$$'})
title({'E2) Radial Quantities'}, 'Interpreter','latex', 'fontsize', fff)
legend('Controls', 'Medulloblastomas', 'fontsize', 8)
hAxes.FontSize = fff2;

% Stats
% --- RK Metric ---
[h, ptotal, ci, stats] = ttest2(rk(controls(:)), rk(tumor(:)));
cohens_d_rk = (mean(rk(controls(:))) - mean(rk(tumor(:)))) / stats.sd;
disp(['RK: ', num2str(cohens_d_rk)]);

% --- VRK Metric ---
[h, pvar, ci, stats] = ttest2(vrk(controls(:)), vrk(tumor(:)));
cohens_d_vrk = (mean(vrk(controls(:))) - mean(vrk(tumor(:)))) / stats.sd;
disp(['VRK: ', num2str(cohens_d_vrk)]);

% --- uRK Metric ---
[h, pmicro, ci, stats] = ttest2(urk(controls(:)), urk(tumor(:)));
cohens_d_urk = (mean(urk(controls(:))) - mean(urk(tumor(:)))) / stats.sd;
disp(['uRK: ', num2str(cohens_d_urk)]);

disp(ptotal)
disp(pvar)
disp(pmicro)

adj_p = [ptotal, pvar, pmicro]*3
raw_pvals = [ptotal, pvar, pmicro];

[h_fdr, crit_p, adj_ci_cvrg, adj_p] = fdr_bh(raw_pvals, 0.05, 'pdep', 'no');

for mi = 1:3
    if mi == 1
        y_max = max(rk(:)) + 0.1; % Adjust 0.1 as needed for spacing
    elseif mi == 2
        y_max = max(vrk(:)) + 0.1;
    elseif mi == 3
        y_max = max(urk(:)) + 0.1;
    end
    
    if adj_p(mi) < 0.001
        sig_marker = '***';
    elseif adj_p(mi) < 0.01
        sig_marker = '**';
    elseif adj_p(mi) < 0.05
        sig_marker = '*';
    else
        sig_marker = '';
    end
    
    if ~isempty(sig_marker)
        text(mi, y_max, sig_marker, ...
            'HorizontalAlignment', 'center', ...
            'VerticalAlignment', 'bottom', ...
            'FontSize', fff, ...
            'Interpreter', 'latex', ...
            'FontWeight', 'bold')
    end
end


%% F3) plot axial kurtosis quatities
hAxes = subplot_tight(3, 4, 11);
hAxes.Position = P3;

ct_all = [ak(controls(:)); vak(controls(:)); uak(controls(:))];
st_all = [ak(tumor(:)); vak(tumor(:)); uak(tumor(:))];
ct = mean(ct_all, 2);
st = mean(st_all, 2);

hAxes.TickLabelInterpreter = 'latex';
b1 = bar(bp-xs, ct, c*x);
b1.FaceColor = ct_color;
hold on
b2 = bar(bp+xs, st, c*x);
b2.FaceColor = st_color;
plot([bp(:)-xs, bp(:)+xs]', [ct, st]', 'black')
plot(sub_control(:)-xs, ct_all(:), '.black')
plot(sub_tumor(:)+xs, st_all(:), '.black')
xlim([0.5, 3.5])
ylim([0, 2])
set(gca,'XTick', [1 2 3], 'YTick',[0 0.5 1 1.5 2 2.5 3]);
xticks([1 2 3])
xticklabels({'$$K^\parallel_t$$', ...
    '$$K^\parallel_v$$',...
    '$$K^\parallel_{\mu}$$'})
title('E3) Axial Quantities', 'Interpreter', 'latex', ...
    'fontsize', fff2)
legend('Controls', 'Medulloblastomas', 'fontsize', 8)
hAxes.FontSize = fff2;

% Stats
% --- AK Metric ---
[h, ptotal, ci, stats] = ttest2(ak(controls(:)), ak(tumor(:)));
cohens_d_ak = (mean(ak(controls(:))) - mean(ak(tumor(:)))) / stats.sd;
disp(['AK: ', num2str(cohens_d_ak)]);

% --- VAK Metric ---
[h, pvar, ci, stats] = ttest2(vak(controls(:)), vak(tumor(:)));
cohens_d_vak = (mean(vak(controls(:))) - mean(vak(tumor(:)))) / stats.sd;
disp(['VAK: ', num2str(cohens_d_vak)]);

% --- uAK Metric ---
[h, pmicro, ci, stats] = ttest2(uak(controls(:)), uak(tumor(:)));
cohens_d_uak = (mean(uak(controls(:))) - mean(uak(tumor(:)))) / stats.sd;
disp(['uAK: ', num2str(cohens_d_uak)]);
disp(ptotal)
disp(pvar)
disp(pmicro)

adj_p = [ptotal, pvar, pmicro]*3
raw_pvals = [ptotal, pvar, pmicro];

[h_fdr, crit_p, adj_ci_cvrg, adj_p] = fdr_bh(raw_pvals, 0.05, 'pdep', 'no');

for mi = 1:3
    if mi == 1
        y_max = max(ak(:)) + 0.1; % Adjust 0.1 as needed for spacing
    elseif mi == 2
        y_max = max(vak(:)) + 0.1;
    elseif mi == 3
        y_max = max(uak(:)) + 0.1;
    end
    
    if adj_p(mi) < 0.001
        sig_marker = '***';
    elseif adj_p(mi) < 0.01
        sig_marker = '**';
    elseif adj_p(mi) < 0.05
        sig_marker = '*';
    else
        sig_marker = '';
    end
    
    if ~isempty(sig_marker)
        text(mi, y_max, sig_marker, ...
            'HorizontalAlignment', 'center', ...
            'VerticalAlignment', 'bottom', ...
            'FontSize', fff, ...
            'Interpreter', 'latex', ...
            'FontWeight', 'bold')
    end
end


%% F4) Anisotropy/dispersion
hAxes = subplot_tight(3, 4, 12);
hAxes.Position = P4;

ct_all = [fa(controls(:)); ufa(controls(:)); op(controls(:))];
st_all = [fa(tumor(:)); ufa(tumor(:)); op(tumor(:))];
ct = mean(ct_all, 2);
st = mean(st_all, 2);

hAxes.TickLabelInterpreter = 'latex';
b1 = bar(bp-xs, ct, c*x);
b1.FaceColor = ct_color;
hold on
b2 = bar(bp+xs, st, c*x);
b2.FaceColor = st_color;
plot([bp(:)-xs, bp(:)+xs]', [ct, st]', 'black')
plot(sub_control(:)-xs, ct_all(:), '.black')
plot(sub_tumor(:)+xs, st_all(:), '.black')
xlim([0.5, 3.5])
ylim([0, 1])
set(gca,'XTick', [1 2 3], 'YTick',[0 0.2 0.4 0.6 0.8 1]);
xticks([1 2 3])
xticklabels({'$$FA$$', ...
    '$${\mu}FA$$',...
    '$$OP$$'})
title('E4) Anisotropy/Dispersion', 'Interpreter', 'latex', ...
    'fontsize', fff2)
legend('Controls', 'Medulloblastomas', 'fontsize', 8)
hAxes.FontSize = fff2;

% Stats
% --- FA Metric ---
[h, ptotal, ci, stats] = ttest2(fa(controls(:)), fa(tumor(:)));
cohens_d_fa = (mean(fa(controls(:))) - mean(fa(tumor(:)))) / stats.sd;
disp(['FA: ', num2str(cohens_d_fa)]);

% --- uFA Metric ---
[h, pvar, ci, stats] = ttest2(ufa(controls(:)), ufa(tumor(:)));
cohens_d_ufa = (mean(ufa(controls(:))) - mean(ufa(tumor(:)))) / stats.sd;
disp(['uFA: ', num2str(cohens_d_ufa)]);

% --- OP Metric ---
[h, pmicro, ci, stats] = ttest2(op(controls(:)), op(tumor(:)));
cohens_d_op = (mean(op(controls(:))) - mean(op(tumor(:)))) / stats.sd;
disp(['OP: ', num2str(cohens_d_op)]);

disp(ptotal)
disp(pvar)
disp(pmicro)

adj_p = [ptotal, pvar, pmicro]*3
raw_pvals = [ptotal, pvar, pmicro];

[h_fdr, crit_p, adj_ci_cvrg, adj_p] = fdr_bh(raw_pvals, 0.05, 'pdep', 'no');

for mi = 1:3
    if mi == 1
        y_max = max(fa(:)) + 0.05;
    elseif mi == 2
        y_max = max(ufa(:)) + 0.05;
    elseif mi == 3
        y_max = max(op(:)) + 0.05;
    end
    
    if adj_p(mi) < 0.001
        sig_marker = '***';
    elseif adj_p(mi) < 0.01
        sig_marker = '**';
    elseif adj_p(mi) < 0.05
        sig_marker = '*';
    else
        sig_marker = '';
    end
    
    if ~isempty(sig_marker)
        text(mi, y_max, sig_marker, ...
            'HorizontalAlignment', 'center', ...
            'VerticalAlignment', 'bottom', ...
            'FontSize', fff, ...
            'Interpreter', 'latex', ...
            'FontWeight', 'bold')
    end
end


% Overarching Title E
annotation('textbox', [0.05, 0.29, 0.90, 0.05], ...
            'String', '$$\textbf{E) ROI Analysis}$$', ...
            'HorizontalAlignment', 'center', ...
            'BackgroundColor', 'none', ...
            'EdgeColor', 'none', 'Interpreter', 'latex', 'fontsize', 12);

print(fig, 'figure5.pdf', '-dpdf', '-painters');
%print(fig, 'figure5.eps', '-depsc', '-painters');
fprintf('Figures saved successfully!\n');




%% Supplementary Figure S2

fig = figure('color', [1 1 1], 'Units', 'centimeters', ...
    'Position', [-30 -5 21 23]);
%figure
set(fig, 'PaperPositionMode', 'auto');
set(fig, 'PaperOrientation', 'portrait');
set(groot, 'DefaultTextFontName', 'Arial');

annotation('line', [0.51 0.51], [0.97 0.35], 'Color', 'k', 'LineWidth', 2)

fff = 10;

ax = subplot_tight(6, 6, 1);
imagesc(squeeze(MDcontrol), [0 2]); axis image; axis off;
hold on
colormap_set
cb.Ticks = [0 2];
text(70, 80, {'A1)  $$\overline{D}$$', '($${\mu}m^2/ms)$$'}, 'BackGroundColor', 'none', 'Interpreter','latex', 'fontsize', fff, 'color', [1 1 1], 'HorizontalAlignment', 'right')


ax = subplot_tight(6, 6, 2);
imagesc(squeeze(RDcontrol), [0 2]); axis image; axis off;
hold on
colormap_set
cb.Ticks = [0 2];
title({'\textbf{A) CTI maps (Control)}'}, 'Interpreter','latex', 'fontsize', fff)
text(70, 80, {'A2)  $$D^\bot$$', '($${\mu}m^2/ms)$$'}, 'BackGroundColor', 'none', 'Interpreter','latex', 'fontsize', fff, 'color', [1 1 1], 'HorizontalAlignment', 'right')

ax = subplot_tight(6, 6, 3);
imagesc(squeeze(ADcontrol), [0 2]); axis image; axis off;
hold on
colormap_set
cb.Ticks = [0 2];
text(70, 80, {'A3)  $$D^\parallel $$', '($${\mu}m^2/ms)$$'}, 'BackGroundColor', 'none', 'Interpreter','latex', 'fontsize', fff, 'color', [1 1 1], 'HorizontalAlignment', 'right')


ax = subplot_tight(6, 6, 4);
imagesc(squeeze(MDtumor), [0 2]); axis image; axis off;
hold on
colormap_set
cb.Ticks = [0 2];
text(70, 85, {'B1)  $$\overline{D}$$'}, 'BackGroundColor', 'none', 'Interpreter','latex', 'fontsize', fff, 'color', [1 1 1], 'HorizontalAlignment', 'right')


ax = subplot_tight(6, 6, 5);
imagesc(squeeze(RDtumor), [0 2]); axis image; axis off;
hold on
colormap_set
cb.Ticks = [0 2];
title({'\textbf{B) CTI maps (Medulloblastoma)}'}, 'Interpreter','latex', 'fontsize', fff)
text(70, 85, {'B2)  $$D^\bot$$'}, 'BackGroundColor', 'none', 'Interpreter','latex', 'fontsize', fff, 'color', [1 1 1], 'HorizontalAlignment', 'right')

ax = subplot_tight(6, 6, 6);
imagesc(squeeze(ADtumor), [0 2]); axis image; axis off;
hold on
colormap_set
cb.Ticks = [0 2];
text(70, 85, {'B3)  $$D^\parallel $$'}, 'BackGroundColor', 'none', 'Interpreter','latex', 'fontsize', fff, 'color', [1 1 1], 'HorizontalAlignment', 'right')


ax = subplot_tight(6, 6, 7);
imagesc(squeeze(MKcontrol), [0  2]); axis image; axis off;
hold on
colormap_set
cb.Ticks = [0 2];
text(70, 85, {'A4)  $$\overline{K}_t$$'}, 'BackGroundColor', 'none', 'Interpreter','latex', 'fontsize', fff, 'color', [1 1 1], 'HorizontalAlignment', 'right')


ax = subplot_tight(6, 6, 8);
imagesc(squeeze(RKcontrol), [0 2]); axis image; axis off;
hold on
colormap_set
cb.Ticks = [0 2];
text(70, 85, {'A5)  $$K_t^\bot$$'}, 'BackGroundColor', 'none', 'Interpreter','latex', 'fontsize', fff, 'color', [1 1 1], 'HorizontalAlignment', 'right')

ax = subplot_tight(6, 6, 9);
imagesc(squeeze(AKcontrol), [0 2]); axis image; axis off;
hold on
colormap_set
cb.Ticks = [0 2];
text(70, 85, {'A6)  $$K_t^\parallel $$'}, 'BackGroundColor', 'none', 'Interpreter','latex', 'fontsize', fff, 'color', [1 1 1], 'HorizontalAlignment', 'right')


ax = subplot_tight(6, 6, 10);
imagesc(squeeze(MKtumor), [0  2]); axis image; axis off;
hold on
colormap_set
cb.Ticks = [0 2];
text(70, 85, {'B4)  $$\overline{K}_t$$'}, 'BackGroundColor', 'none', 'Interpreter','latex', 'fontsize', fff, 'color', [1 1 1], 'HorizontalAlignment', 'right')


ax = subplot_tight(6, 6, 11);
imagesc(squeeze(RKtumor), [0 2]); axis image; axis off;
hold on
colormap_set
cb.Ticks = [0 2];
text(70, 85, {'B5)  $$K_t^\bot$$'}, 'BackGroundColor', 'none', 'Interpreter','latex', 'fontsize', fff, 'color', [1 1 1], 'HorizontalAlignment', 'right')

ax = subplot_tight(6, 6, 12);
imagesc(squeeze(AKtumor), [0 2]); axis image; axis off;
hold on
colormap_set
cb.Ticks = [0 2];
text(70, 85, {'B6)  $$K_t^\parallel $$'}, 'BackGroundColor', 'none', 'Interpreter','latex', 'fontsize', fff, 'color', [1 1 1], 'HorizontalAlignment', 'right')



ax = subplot_tight(6, 6, 13);
imagesc(squeeze(MKvcontrol), [0  1.5]); axis image; axis off;
hold on
colormap_set
cb.Ticks = [0 1.5];
text(70, 85, {'A7)  $$\overline{K}_v$$'}, 'BackGroundColor', 'none', 'Interpreter','latex', 'fontsize', fff, 'color', [1 1 1], 'HorizontalAlignment', 'right')


ax = subplot_tight(6, 6, 14);
imagesc(squeeze(RKvcontrol), [0 1.5]); axis image; axis off;
hold on
colormap_set
cb.Ticks = [0 1.5];
text(70, 85, {'A8)  $$K_v^\bot$$'}, 'BackGroundColor', 'none', 'Interpreter','latex', 'fontsize', fff, 'color', [1 1 1], 'HorizontalAlignment', 'right')

ax = subplot_tight(6, 6, 15);
imagesc(squeeze(AKvcontrol), [0 1.5]); axis image; axis off;
hold on
colormap_set
cb.Ticks = [0 1.5];
text(70, 85, {'A9)  $$K_v^\parallel $$'}, 'BackGroundColor', 'none', 'Interpreter','latex', 'fontsize', fff, 'color', [1 1 1], 'HorizontalAlignment', 'right')


ax = subplot_tight(6, 6, 16);
imagesc(squeeze(MKvtumor), [0  2]); axis image; axis off;
hold on
colormap_set
cb.Ticks = [0 2];
text(70, 85, {'B7)  $$\overline{K}_v$$'}, 'BackGroundColor', 'none', 'Interpreter','latex', 'fontsize', fff, 'color', [1 1 1], 'HorizontalAlignment', 'right')


ax = subplot_tight(6, 6, 17);
imagesc(squeeze(RKvtumor), [0 2]); axis image; axis off;
hold on
colormap_set
cb.Ticks = [0 2];
text(70, 85, {'B8)  $$K_v^\bot$$'}, 'BackGroundColor', 'none', 'Interpreter','latex', 'fontsize', fff, 'color', [1 1 1], 'HorizontalAlignment', 'right')

ax = subplot_tight(6, 6, 18);
imagesc(squeeze(AKvtumor), [0 2]); axis image; axis off;
hold on
colormap_set
cb.Ticks = [0 2];
text(70, 85, {'B9)  $$K_v^\parallel $$'}, 'BackGroundColor', 'none', 'Interpreter','latex', 'fontsize', fff, 'color', [1 1 1], 'HorizontalAlignment', 'right')


ax = subplot_tight(6, 6, 20);
imagesc(squeeze(OPcontrol), [0 1]); axis image; axis off;
ax.Position(1) = ax.Position(1) + 0.08; 
hold on
colormap_set
cb.Ticks = [0 1];
text(70, 85, {'A11)  $$OP$$'}, 'BackGroundColor', 'none', 'Interpreter','latex', 'fontsize', fff, 'color', [1 1 1], 'HorizontalAlignment', 'right')


ax = subplot_tight(6, 6, 19);
imagesc(squeeze(uFAcontrol), [0 1]); axis image; axis off;
ax.Position(1) = ax.Position(1) + 0.08; 
hold on
colormap_set
cb.Ticks = [0 1];
text(70, 85, {'A10)  $${\mu}FA$$'}, 'BackGroundColor', 'none', 'Interpreter','latex', 'fontsize', fff, 'color', [1 1 1], 'HorizontalAlignment', 'right')


ax = subplot_tight(6, 6, 23);
imagesc(squeeze(OPtumor), [0 1]); axis image; axis off;
ax.Position(1) = ax.Position(1) + 0.08; 
hold on
colormap_set
cb.Ticks = [0 1];
text(70, 85, {'B11)  $$OP$$'}, 'BackGroundColor', 'none', 'Interpreter','latex', 'fontsize', fff, 'color', [1 1 1], 'HorizontalAlignment', 'right')


ax = subplot_tight(6, 6, 22);
imagesc(squeeze(uFAtumor), [0 1]); axis image; axis off;
ax.Position(1) = ax.Position(1) + 0.08; 
hold on
colormap_set
cb.Ticks = [0 1];
text(70, 85, {'B10)  $${\mu}FA$$'}, 'BackGroundColor', 'none', 'Interpreter','latex', 'fontsize', fff, 'color', [1 1 1], 'HorizontalAlignment', 'right')



hAxes = subplot(4, 3, 11);
hAxes.Position = [0.35 0.1 0.3 0.2];

ct_all = [md(controls(:)); rd(controls(:)); ad(controls(:))];
st_all = [md(tumor(:)); rd(tumor(:)); ad(tumor(:))];
ct = mean(ct_all, 2);
st = mean(st_all, 2);

hAxes.TickLabelInterpreter = 'latex';
b1 = bar(bp-xs, ct, c*x);
b1.FaceColor = ct_color;
hold on
b2 = bar(bp+xs, st, c*x);
b2.FaceColor = st_color;
plot([bp(:)-xs, bp(:)+xs]', [ct, st]', 'black')
plot(sub_control(:)-xs, ct_all(:), '.black')
plot(sub_tumor(:)+xs, st_all(:), '.black')
xlim([0.5, 3.5])
ylim([0, 1])
set(gca,'XTick', [1 2 3], 'YTick',[0 0.2 0.4 0.6 0.8 1]);
xticks([1 2 3])

xticklabels({'$$\overline{D}$$', ...
    '$$D^\bot_v$$',...
    '$$D^\parallel$$'})

title('C) Diffusion Quantities', 'fontsize', fff2)
legend('Controls', 'Medulloblastomas', 'fontsize', 8)
hAxes.FontSize = fff2;

% Stats
% --- MD Metric ---
[h, ptotal, ci, stats] = ttest2(md(controls(:)), md(tumor(:)));
cohens_d_md = (mean(md(controls(:))) - mean(md(tumor(:)))) / stats.sd;
disp(['MD: ', num2str(cohens_d_md)]);

% --- RD Metric ---
[h, pvar, ci, stats] = ttest2(rd(controls(:)), rd(tumor(:)));
cohens_d_rd = (mean(rd(controls(:))) - mean(rd(tumor(:)))) / stats.sd;
disp(['RD: ', num2str(cohens_d_rd)]);

% --- AD Metric ---
[h, pmicro, ci, stats] = ttest2(ad(controls(:)), ad(tumor(:)));
cohens_d_ad = (mean(ad(controls(:))) - mean(ad(tumor(:)))) / stats.sd;
disp(['AD: ', num2str(cohens_d_ad)]);
disp(ptotal)
disp(pvar)
disp(pmicro)



adj_p = [ptotal, pvar, pmicro]*3
raw_pvals = [ptotal, pvar, pmicro];

[h_fdr, crit_p, adj_ci_cvrg, adj_p] = fdr_bh(raw_pvals, 0.05, 'pdep', 'no');

for mi = 1:3
    if mi == 1
        y_max = max(md(:)) + 0.05;
    elseif mi == 2
        y_max = max(rd(:)) + 0.05;
    elseif mi == 3
        y_max = max(ad(:)) + 0.05;
    end
    
    if adj_p(mi) < 0.001
        sig_marker = '***';
    elseif adj_p(mi) < 0.01
        sig_marker = '**';
    elseif adj_p(mi) < 0.05
        sig_marker = '*';
    else
        sig_marker = '';
    end
    
    if ~isempty(sig_marker)
        text(mi, y_max, sig_marker, ...
            'HorizontalAlignment', 'center', ...
            'VerticalAlignment', 'bottom', ...
            'FontSize', fff, ...
            'Interpreter', 'latex', ...
            'FontWeight', 'bold')
    end
end

print(fig, 'figureS2.pdf', '-dpdf', '-painters');
%print(fig, 'figure5.eps', '-depsc', '-painters');
fprintf('Figures saved successfully!\n');

