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
annotation('textbox', [0.05, 0.82, 0.4, 0.1], ...
            'String', '$$\textbf{A) Histology (CT2A)}$$', ...
            'HorizontalAlignment', 'center', ...
            'BackgroundColor', 'none', ...
            'EdgeColor', 'none', 'Interpreter', 'latex', 'fontsize', fff);

annotation('textbox', [0.52, 0.82, 0.4, 0.1], ...
            'String', '$$\textbf{B) Histology (GL261)}$$', ...
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
axis equal
axis off
text(122, 50, 'A1)', 'BackGroundColor', [0 0 0], 'Interpreter','latex', 'fontsize', 10, 'color', [1 1 1], 'HorizontalAlignment', 'center')
ax1.Position = P1;

ax1 = subplot_tight(3, 4, 2);
Image2 = imread('panelA2.jpg');
imagesc(Image2)
axis equal
axis off
text(122, 50, 'A2)', 'BackGroundColor', [0 0 0], 'Interpreter','latex', 'fontsize', 10, 'color', [1 1 1], 'HorizontalAlignment', 'center')
ax1.Position = P2;

ax1 = subplot_tight(3, 4, 3);
Image2 = imread('panelB1.jpg');
imagesc(Image2)
axis equal
axis off
text(122, 50, 'B1)', 'BackGroundColor', [0 0 0], 'Interpreter','latex', 'fontsize', 10, 'color', [1 1 1], 'HorizontalAlignment', 'center')
ax1.Position = P3;

ax1 = subplot_tight(3, 4, 4);
Image2 = imread('panelB2.jpg');
imagesc(Image2)
axis equal
axis off
text(122, 50, 'B2)', 'BackGroundColor', [0 0 0], 'Interpreter','latex', 'fontsize', 10, 'color', [1 1 1], 'HorizontalAlignment', 'center')
ax1.Position = P4;

%% Analysis 
rnh_005_paper_analysis
fff = 12;
fff2 = 10;

%% Panels C
mycmap = turbo(1000);

%FA
ax = subplot_tight(6, 6, 13);
imagesc(squeeze(FAct2a), [0 1]); axis image; axis off;
hold on
%title('C1) $$FA$$', 'Interpreter','latex', 'fontsize', fff)
colormap_set
cb.Ticks = [0 1];
text(30, 8, 'C1)  $$FA$$', 'BackGroundColor', 'none', 'Interpreter','latex', 'fontsize', fff2, 'color', [1 1 1], 'HorizontalAlignment', 'center')

% KANISO
ax = subplot_tight(6, 6, 14);
imagesc(squeeze(KANIct2a), [0 1.5]); axis image; axis off;
hold on
%title({'\textbf{C) CTI maps (CT2A)}'; 'C2) $$\overline{K}_{ani}$$'}, 'Interpreter','latex', 'fontsize', fff)
colormap_set
cb.Ticks = [0 1.5];
text(30, 8, 'C2) $$\overline{K}_{ani}$$', 'BackGroundColor', 'none', 'Interpreter','latex', 'fontsize', fff2, 'color', [1 1 1], 'HorizontalAlignment', 'center')
%title({'\textbf{C) CTI maps (CT2A)}'}, 'Interpreter','latex', 'fontsize', fff, 'HorizontalAlignment', 'l'),...


% Kiso
ax = subplot_tight(6, 6, 15);
imagesc(squeeze(KISOct2a), [0 1.5]); axis image; axis off;
hold on
%title({'C3) $$\overline{K}_{iso}$$'}, 'Interpreter','latex', 'fontsize', fff)
colormap_set
cb.Ticks = [0 1.5];
text(30, 8, 'C3) $$\overline{K}_{iso}$$', 'BackGroundColor', 'none', 'Interpreter','latex', 'fontsize', fff2, 'color', [1 1 1], 'HorizontalAlignment', 'center')

% Ku
ax = subplot_tight(6, 6, 16);
imagesc(squeeze(MKict2a), [0 1.5]); axis image; axis off;
hold on
%title('C4) $$\overline{K}_{\mu}$$', 'Interpreter','latex', 'fontsize', fff)
colormap_set
cb.Ticks = [0 1.5];
text(30, 8, 'C4) $$\overline{K}_{\mu}$$', 'BackGroundColor', 'none', 'Interpreter','latex', 'fontsize', fff2, 'color', [1 1 1], 'HorizontalAlignment', 'center')
%title({'\textbf{C) CTI maps (CT2A)}'}, 'Interpreter','latex', 'fontsize', fff, 'HorizontalAlignment', 'r')


% RKi
ax = subplot_tight(6, 6, 17);
imagesc(squeeze(RKict2a), [0 1.5]); axis image; axis off;
hold on
%title('C5) $$K^\bot_{\mu}$$', 'Interpreter','latex', 'fontsize', fff)
colormap_set
cb.Ticks = [0 1.5];
text(30, 8, 'C5) $$K^\bot_{\mu}$$', 'BackGroundColor', 'none', 'Interpreter','latex', 'fontsize', fff2, 'color', [1 1 1], 'HorizontalAlignment', 'center')


% AKi
ax = subplot_tight(6, 6, 18);
imagesc(squeeze(AKict2a), [0 1.5]); axis image; axis off;
hold on
%title('C6) $$K^\parallel_{\mu}$$', 'Interpreter','latex', 'fontsize', fff)
colormap_set
cb.Ticks = [0 1.5];
text(30, 8, 'C6) $$K^\parallel_{\mu}$$', 'BackGroundColor', 'none', 'Interpreter','latex', 'fontsize', fff2, 'color', [1 1 1], 'HorizontalAlignment', 'center')

% Overarching Title C
annotation('textbox', [0.05, 0.62, 0.90, 0.05], ...
            'String', '$$\textbf{C) CTI maps (CT2A)}$$', ...
            'HorizontalAlignment', 'center', ...
            'BackgroundColor', 'none', ...
            'EdgeColor', 'none', 'Interpreter', 'latex', 'fontsize', fff);

%% Panels D
%FA
ax = subplot_tight(6, 6, 19);
imagesc(squeeze(FAgl261), [0 1]); axis image; axis off;
hold on
%title('D1) $$FA$$', 'Interpreter','latex', 'fontsize', fff)
colormap_set
cb.Ticks = [0 1];
text(30, 8, 'D1) $$FA$$', 'BackGroundColor', 'none', 'Interpreter','latex', 'fontsize', fff2, 'color', [1 1 1], 'HorizontalAlignment', 'center')

% KANISO
ax = subplot_tight(6, 6, 20);
imagesc(squeeze(KANIgl261), [0 1.5]); axis image; axis off;
hold on
%title({'\textbf{D) CTI maps (GL261)}'; 'D2) $$\overline{K}_{ani}$$'}, 'Interpreter','latex', 'fontsize', fff)
colormap_set
cb.Ticks = [0 1.5];
text(30, 8, 'D2) $$\overline{K}_{ani}$$', 'BackGroundColor', 'none', 'Interpreter','latex', 'fontsize', fff2, 'color', [1 1 1], 'HorizontalAlignment', 'center')

% Kiso
ax = subplot_tight(6, 6, 21);
imagesc(squeeze(KISOgl261), [0 1.5]); axis image; axis off;
hold on
%title({'D3) $$\overline{K}_{iso}$$'}, 'Interpreter','latex', 'fontsize', fff)
colormap_set
cb.Ticks = [0 1.5];
text(30, 8, 'D3) $$\overline{K}_{iso}$$', 'BackGroundColor', 'none', 'Interpreter','latex', 'fontsize', fff2, 'color', [1 1 1], 'HorizontalAlignment', 'center')

% Ku
ax = subplot_tight(6, 6, 22);
imagesc(squeeze(MKigl261), [0 1.5]); axis image; axis off;
hold on
%title('D4) $$\overline{K}_{\mu}$$', 'Interpreter','latex', 'fontsize', fff)
colormap_set
cb.Ticks = [0 1.5];
text(30, 8, 'D4) $$\overline{K}_{\mu}$$', 'BackGroundColor', 'none', 'Interpreter','latex', 'fontsize', fff2, 'color', [1 1 1], 'HorizontalAlignment', 'center')

%title({'\textbf{D) CTI maps (GL261)}'}, 'Interpreter','latex', 'fontsize', fff, 'HorizontalAlignment', 'r')

% RKi
ax = subplot_tight(6, 6, 23);
imagesc(squeeze(RKigl261), [0 1.5]); axis image; axis off;
hold on
%title('D5) $$K^\bot_{\mu}$$', 'Interpreter','latex', 'fontsize', fff)
colormap_set
cb.Ticks = [0 1.5];
text(30, 8, 'D5) $$K^\bot_{\mu}$$', 'BackGroundColor', 'none', 'Interpreter','latex', 'fontsize', fff2, 'color', [1 1 1], 'HorizontalAlignment', 'center')

% AKi
ax = subplot_tight(6, 6, 24);
imagesc(squeeze(AKigl261), [0 1.5]); axis image; axis off;
hold on
%title('D6) $$K^\parallel_{\mu}$$', 'Interpreter','latex', 'fontsize', fff)
colormap_set
cb.Ticks = [0 1.5];
text(30, 8, 'D6) $$K^\parallel_{\mu}$$', 'BackGroundColor', 'none', 'Interpreter','latex', 'fontsize', fff2, 'color', [1 1 1], 'HorizontalAlignment', 'center')


% Overarching Title D
annotation('textbox', [0.05, 0.46, 0.90, 0.05], ...
            'String', '$$\textbf{D) CTI maps (GL261)}$$', ...
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
fff = 14;

ct_color = [0.5 0.5 0.5];
st_color = [1 1 1];

bp = 1:4;
sub_control = [ones(1, 3); 2*ones(1, 3); 3*ones(1, 3); 4*ones(1, 3)];
sub_tumor = [ones(1, 3); 2*ones(1, 3); 3*ones(1, 3); 4*ones(1, 3)];

%% E1) plot mean kurtosis quatities
hAxes = subplot(3, 4, 9);
hAxes.Position = P1;

mk = [];
kaniso = [];
kiso = [];
uk = [];
for ti=1:6
    mk = [mk, mean(gliomas(ti).MKv)];
    kaniso = [kaniso, mean(gliomas(ti).KANISOv)];
    kiso = [kiso, mean(gliomas(ti).KISOv)];
    uk = [kiso, mean(gliomas(ti).KISOv)];
end
ct_all = [mk(1:3); kaniso(1:3); kiso(1:3); uk(1:3)];
st_all = [mk(4:6); kaniso(4:6); kiso(4:6); uk(4:6)];
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
ylim([0, 1])
set(gca,'XTick', [1 2 3 4], 'YTick',[0 0.5 1 1.5 2 2.5 3]);
xticks([1 2 3 4])
xticklabels({'$$\overline{K}_t$$', ...
    '$$\overline{K}_{ani}$$',...
    '$$\overline{K}_{iso}$$',...
    '$$\overline{K}_{\mu}$$'})
title('E1) Mean Quantities', 'Interpreter', 'latex', ...
    'fontsize', fff2)
legend('CT2A', 'GL261', 'fontsize', fff2)
hAxes.FontSize = fff2;

% Stats
[h,ptotal,ci,stats] = ttest2(mk(1:3), mk(4:6));
cohens_d_mk = (mean(mk(1:3)) - mean(mk(4:6))) / stats.sd;
disp(['MK Cohen''s d: ', num2str(cohens_d_mk)]);

[h,pvar1,ci,stats] =  ttest2(kaniso(1:3), kaniso(4:6));
cohens_d_kaniso = (mean(kaniso(1:3)) - mean(kaniso(4:6))) / stats.sd;
disp(['K_aniso Cohen''s d: ', num2str(cohens_d_kaniso)]);

[h,pvar2,ci,stats] = ttest2(kiso(1:3), kiso(4:6));
cohens_d_kiso = (mean(kiso(1:3)) - mean(kiso(4:6))) / stats.sd;
disp(['K_iso Cohen''s d: ', num2str(cohens_d_kiso)]);

[h,pmicro,ci,stats] = ttest2(uk(1:3), uk(4:6));
cohens_d_uk = (mean(uk(1:3)) - mean(uk(4:6))) / stats.sd;
disp(['uK Cohen''s d: ', num2str(cohens_d_uk)]);
disp(ptotal)
disp(pvar1)
disp(pvar2)
disp(pmicro)

adj_p = [ptotal, pvar1, pvar2, pmicro];

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
sub_control = [ones(1, 3); 2*ones(1, 3); 3*ones(1, 3)];
sub_tumor = [ones(1, 3); 2*ones(1, 3); 3*ones(1, 3)];

hAxes = subplot_tight(3, 4, 10);
hAxes.Position = P2;

rk = [];
vrk = [];
urk = [];
for ti=1:6
    rk = [rk, mean(gliomas(ti).RKv)];
    vrk = [vrk, mean(gliomas(ti).RKTvv)];
    urk = [urk, mean(gliomas(ti).RKTiv)];
end
ct_all = [rk(1:3); vrk(1:3); urk(1:3)];
st_all = [rk(4:6); vrk(4:6); urk(4:6)];
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
set(gca,'XTick', [1 2 3], 'YTick',[0 0.5 1 1.5 2 2.5 3]);
xticks([1 2 3])
xticklabels({'$$K^\bot_t$$', ...
    '$$K^\bot_v$$',...
    '$$K^\bot_{\mu}$$'})
title({'E2) Radial Quantities'}, 'Interpreter','latex', 'fontsize', fff)
legend('CT2A', 'GL261', 'fontsize', fff2)
hAxes.FontSize = fff2;

% Stats
[h,ptotal,ci,stats] = ttest2(rk(1:3), rk(4:6));
cohens_d_rk = (mean(rk(1:3)) - mean(rk(4:6))) / stats.sd;
disp(['RK Cohen''s d: ', num2str(cohens_d_rk)]);

[h,pvar,ci,stats] = ttest2(vrk(1:3), vrk(4:6));
cohens_d_vrk = (mean(vrk(1:3)) - mean(vrk(4:6))) / stats.sd;
disp(['VRK Cohen''s d: ', num2str(cohens_d_vrk)]);

[h,pmicro,ci,stats] = ttest2(urk(1:3), urk(4:6));
cohens_d_urk = (mean(urk(1:3)) - mean(urk(4:6))) / stats.sd;
disp(['uRK Cohen''s d: ', num2str(cohens_d_urk)]);
disp(ptotal)
disp(pvar)
disp(pmicro)

adj_p = [ptotal, pvar, pmicro];

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

ak = [];
vak = [];
uak = [];
for ti=1:6
    ak = [ak, mean(gliomas(ti).AKv)];
    vak = [vak, mean(gliomas(ti).AKTvv)];
    uak = [uak, mean(gliomas(ti).AKTiv)];
end
ct_all = [ak(1:3); vak(1:3); uak(1:3)];
st_all = [ak(4:6); vak(4:6); uak(4:6)];
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
set(gca,'XTick', [1 2 3], 'YTick',[0 0.5 1 1.5 2 2.5 3]);
xticks([1 2 3])
xticklabels({'$$K^\parallel_t$$', ...
    '$$K^\parallel_v$$',...
    '$$K^\parallel_{\mu}$$'})
title('E3) Axial Quantities', 'Interpreter', 'latex', ...
    'fontsize', fff2)
legend('CT2A', 'GL261', 'fontsize', fff2)
hAxes.FontSize = fff2;

% Stats
[h,ptotal,ci,stats] = ttest2(ak(1:3), ak(4:6));
cohens_d_ak = (mean(ak(1:3)) - mean(ak(4:6))) / stats.sd;
disp(['AK Cohen''s d: ', num2str(cohens_d_ak)]);

[h,pvar,ci,stats] = ttest2(vak(1:3), vak(4:6));
cohens_d_vak = (mean(vak(1:3)) - mean(vak(4:6))) / stats.sd;
disp(['VAK Cohen''s d: ', num2str(cohens_d_vak)]);

[h,pmicro,ci,stats] = ttest2(uak(1:3), uak(4:6));
cohens_d_uak = (mean(uak(1:3)) - mean(uak(4:6))) / stats.sd;
disp(['uAK Cohen''s d: ', num2str(cohens_d_uak)]);
disp(ptotal)
disp(pvar)
disp(pmicro)

adj_p = [ptotal, pvar, pmicro];

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

ak = [];
vak = [];
uak = [];
for ti=1:6
    ak = [ak, mean(gliomas(ti).FAv)];
    vak = [vak, mean(gliomas(ti).uFAv)];
    uak = [uak, mean(gliomas(ti).OPv)];
end
ct_all = [ak(1:3); vak(1:3); uak(1:3)];
st_all = [ak(4:6); vak(4:6); uak(4:6)];
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
legend('CT2A', 'GL261', 'fontsize', fff2, 'Location', 'northwest')
hAxes.FontSize = fff2;

% Stats
[h,ptotal,ci,stats] = ttest2(ak(1:3), ak(4:6));
cohens_d_fa = (mean(ak(1:3)) - mean(ak(4:6))) / stats.sd;
disp(['FA Cohen''s d: ', num2str(cohens_d_fa)]);

[h,pvar,ci,stats] = ttest2(vak(1:3), vak(4:6));
cohens_d_ufa = (mean(vak(1:3)) - mean(vak(4:6))) / stats.sd;
disp(['uFA Cohen''s d: ', num2str(cohens_d_ufa)]);

[h,pmicro,ci,stats] = ttest2(uak(1:3), uak(4:6));
cohens_d_op = (mean(uak(1:3)) - mean(uak(4:6))) / stats.sd;
disp(['OP Cohen''s d: ', num2str(cohens_d_op)]);
disp(ptotal)
disp(pvar)
disp(pmicro)

adj_p = [ptotal, pvar, pmicro]

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

annotation('textbox', [0.05, 0.29, 0.90, 0.05], ...
            'String', '$$\textbf{E) ROI Analysis}$$', ...
            'HorizontalAlignment', 'center', ...
            'BackgroundColor', 'none', ...
            'EdgeColor', 'none', 'Interpreter', 'latex', 'fontsize', 12);


print(fig, 'figure6.pdf', '-dpdf', '-painters');
%print(fig, 'figure6.eps', '-depsc', '-painters');
fprintf('Figures saved successfully!\n');


%% Supplementary Figure S3

fig = figure('color', [1 1 1], 'Units', 'centimeters', ...
    'Position', [-30 -5 21 23]);
%figure
set(fig, 'PaperPositionMode', 'auto');
set(fig, 'PaperOrientation', 'portrait');
set(groot, 'DefaultTextFontName', 'Arial');

annotation('line', [0.51 0.51], [0.97 0.35], 'Color', 'k', 'LineWidth', 2)

fff = 10;

ax = subplot_tight(6, 6, 1);
imagesc(squeeze(MDct2a), [0 2]); axis image; axis off;
hold on
colormap_set
cb.Ticks = [0 2];
text(30, 8, {'A1)  $$\overline{D}$$'}, 'BackGroundColor', 'none', 'Interpreter','latex', 'fontsize', fff, 'color', [1 1 1], 'HorizontalAlignment', 'center')


ax = subplot_tight(6, 6, 2);
imagesc(squeeze(RDct2a), [0 2]); axis image; axis off;
hold on
colormap_set
cb.Ticks = [0 2];
title({'\textbf{A) CTI maps (CT2A)}'}, 'Interpreter','latex', 'fontsize', fff)
text(30, 8, {'A2)  $$D^\bot$$'}, 'BackGroundColor', 'none', 'Interpreter','latex', 'fontsize', fff, 'color', [1 1 1], 'HorizontalAlignment', 'center')

ax = subplot_tight(6, 6, 3);
imagesc(squeeze(ADct2a), [0 2]); axis image; axis off;
hold on
colormap_set
cb.Ticks = [0 2];
text(30, 8, {'A3)  $$D^\parallel $$'}, 'BackGroundColor', 'none', 'Interpreter','latex', 'fontsize', fff, 'color', [1 1 1], 'HorizontalAlignment', 'center')


ax = subplot_tight(6, 6, 4);
imagesc(squeeze(MDgl261), [0 2]); axis image; axis off;
hold on
colormap_set
cb.Ticks = [0 2];
text(30, 8, {'B1)  $$\overline{D}$$'}, 'BackGroundColor', 'none', 'Interpreter','latex', 'fontsize', fff, 'color', [1 1 1], 'HorizontalAlignment', 'center')


ax = subplot_tight(6, 6, 5);
imagesc(squeeze(RDgl261), [0 2]); axis image; axis off;
hold on
colormap_set
cb.Ticks = [0 2];
title({'\textbf{B) CTI maps (Medulloblastoma)}'}, 'Interpreter','latex', 'fontsize', fff)
text(30, 8, {'B2)  $$D^\bot$$'}, 'BackGroundColor', 'none', 'Interpreter','latex', 'fontsize', fff, 'color', [1 1 1], 'HorizontalAlignment', 'center')

ax = subplot_tight(6, 6, 6);
imagesc(squeeze(ADgl261), [0 2]); axis image; axis off;
hold on
colormap_set
cb.Ticks = [0 2];
text(30, 8, {'B3)  $$D^\parallel $$'}, 'BackGroundColor', 'none', 'Interpreter','latex', 'fontsize', fff, 'color', [1 1 1], 'HorizontalAlignment', 'center')


ax = subplot_tight(6, 6, 7);
imagesc(squeeze(MKct2a), [0  2]); axis image; axis off;
hold on
colormap_set
cb.Ticks = [0 2];
text(30, 8, {'A4)  $$\overline{K}_t$$'}, 'BackGroundColor', 'none', 'Interpreter','latex', 'fontsize', fff, 'color', [1 1 1], 'HorizontalAlignment', 'center')


ax = subplot_tight(6, 6, 8);
imagesc(squeeze(RKct2a), [0 2]); axis image; axis off;
hold on
colormap_set
cb.Ticks = [0 2];
text(30, 8, {'A5)  $$K_t^\bot$$'}, 'BackGroundColor', 'none', 'Interpreter','latex', 'fontsize', fff, 'color', [1 1 1], 'HorizontalAlignment', 'center')

ax = subplot_tight(6, 6, 9);
imagesc(squeeze(AKct2a), [0 2]); axis image; axis off;
hold on
colormap_set
cb.Ticks = [0 2];
text(30, 8, {'A6)  $$K_t^\parallel $$'}, 'BackGroundColor', 'none', 'Interpreter','latex', 'fontsize', fff, 'color', [1 1 1], 'HorizontalAlignment', 'center')


ax = subplot_tight(6, 6, 10);
imagesc(squeeze(MKgl261), [0  2]); axis image; axis off;
hold on
colormap_set
cb.Ticks = [0 2];
text(30, 8, {'B4)  $$\overline{K}_t$$'}, 'BackGroundColor', 'none', 'Interpreter','latex', 'fontsize', fff, 'color', [1 1 1], 'HorizontalAlignment', 'center')


ax = subplot_tight(6, 6, 11);
imagesc(squeeze(RKgl261), [0 2]); axis image; axis off;
hold on
colormap_set
cb.Ticks = [0 2];
text(30, 8, {'B5)  $$K_t^\bot$$'}, 'BackGroundColor', 'none', 'Interpreter','latex', 'fontsize', fff, 'color', [1 1 1], 'HorizontalAlignment', 'center')

ax = subplot_tight(6, 6, 12);
imagesc(squeeze(AKgl261), [0 2]); axis image; axis off;
hold on
colormap_set
cb.Ticks = [0 2];
text(30, 8, {'B6)  $$K_t^\parallel $$'}, 'BackGroundColor', 'none', 'Interpreter','latex', 'fontsize', fff, 'color', [1 1 1], 'HorizontalAlignment', 'center')



ax = subplot_tight(6, 6, 13);
imagesc(squeeze(MKvct2a), [0  1.5]); axis image; axis off;
hold on
colormap_set
cb.Ticks = [0 1.5];
text(30, 8, {'A7)  $$\overline{K}_v$$'}, 'BackGroundColor', 'none', 'Interpreter','latex', 'fontsize', fff, 'color', [1 1 1], 'HorizontalAlignment', 'center')


ax = subplot_tight(6, 6, 14);
imagesc(squeeze(RKvct2a), [0 1.5]); axis image; axis off;
hold on
colormap_set
cb.Ticks = [0 1.5];
text(30, 8, {'A8)  $$K_v^\bot$$'}, 'BackGroundColor', 'none', 'Interpreter','latex', 'fontsize', fff, 'color', [1 1 1], 'HorizontalAlignment', 'center')

ax = subplot_tight(6, 6, 15);
imagesc(squeeze(AKvct2a), [0 1.5]); axis image; axis off;
hold on
colormap_set
cb.Ticks = [0 1.5];
text(30, 8, {'A9)  $$K_v^\parallel $$'}, 'BackGroundColor', 'none', 'Interpreter','latex', 'fontsize', fff, 'color', [1 1 1], 'HorizontalAlignment', 'center')


ax = subplot_tight(6, 6, 16);
imagesc(squeeze(MKvgl261), [0  2]); axis image; axis off;
hold on
colormap_set
cb.Ticks = [0 2];
text(30, 8, {'B7)  $$\overline{K}_v$$'}, 'BackGroundColor', 'none', 'Interpreter','latex', 'fontsize', fff, 'color', [1 1 1], 'HorizontalAlignment', 'center')


ax = subplot_tight(6, 6, 17);
imagesc(squeeze(RKvgl261), [0 2]); axis image; axis off;
hold on
colormap_set
cb.Ticks = [0 2];
text(30, 8, {'B8)  $$K_v^\bot$$'}, 'BackGroundColor', 'none', 'Interpreter','latex', 'fontsize', fff, 'color', [1 1 1], 'HorizontalAlignment', 'center')

ax = subplot_tight(6, 6, 18);
imagesc(squeeze(AKvgl261), [0 2]); axis image; axis off;
hold on
colormap_set
cb.Ticks = [0 2];
text(30, 8, {'B9)  $$K_v^\parallel $$'}, 'BackGroundColor', 'none', 'Interpreter','latex', 'fontsize', fff, 'color', [1 1 1], 'HorizontalAlignment', 'center')


ax = subplot_tight(6, 6, 20);
imagesc(squeeze(OPct2a), [0 1]); axis image; axis off;
ax.Position(1) = ax.Position(1) + 0.08; 
hold on
colormap_set
cb.Ticks = [0 1];
text(30, 8, {'A11)  $$OP$$'}, 'BackGroundColor', 'none', 'Interpreter','latex', 'fontsize', fff, 'color', [1 1 1], 'HorizontalAlignment', 'center')


ax = subplot_tight(6, 6, 19);
imagesc(squeeze(uFAct2a), [0 1]); axis image; axis off;
ax.Position(1) = ax.Position(1) + 0.08; 
hold on
colormap_set
cb.Ticks = [0 1];
text(30, 8, {'A10)  $${\mu}FA$$'}, 'BackGroundColor', 'none', 'Interpreter','latex', 'fontsize', fff, 'color', [1 1 1], 'HorizontalAlignment', 'center')


ax = subplot_tight(6, 6, 23);
imagesc(squeeze(OPgl261), [0 1]); axis image; axis off;
ax.Position(1) = ax.Position(1) + 0.08; 
hold on
colormap_set
cb.Ticks = [0 1];
text(30, 8, {'B11)  $$OP$$'}, 'BackGroundColor', 'none', 'Interpreter','latex', 'fontsize', fff, 'color', [1 1 1], 'HorizontalAlignment', 'center')


ax = subplot_tight(6, 6, 22);
imagesc(squeeze(uFAgl261), [0 1]); axis image; axis off;
ax.Position(1) = ax.Position(1) + 0.08; 
hold on
colormap_set
cb.Ticks = [0 1];
text(30, 8, {'B10)  $${\mu}FA$$'}, 'BackGroundColor', 'none', 'Interpreter','latex', 'fontsize', fff, 'color', [1 1 1], 'HorizontalAlignment', 'center')




hAxes = subplot_tight(4, 3, 11);
hAxes.Position = [0.35 0.1 0.3 0.2];

ak = [];
vak = [];
uak = [];
for ti=1:6
    ak = [ak, mean(gliomas(ti).MDv)];
    vak = [vak, mean(gliomas(ti).RDv)];
    uak = [uak, mean(gliomas(ti).ADv)];
end
ct_all = [ak(1:3); vak(1:3); uak(1:3)];
st_all = [ak(4:6); vak(4:6); uak(4:6)];
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
ylim([0, 3])
set(gca,'XTick', [1 2 3], 'YTick',[0 0.5 1 1.5 2 2.5 3]);
xticks([1 2 3])
xticklabels({'$$\overline{D}$$', ...
    '$$D^\bot_v$$',...
    '$$D^\parallel$$'})
title('C) Diffusion Quantities', 'Interpreter', 'latex', ...
    'fontsize', fff2)
legend('CT2A', 'GL261', 'fontsize', fff2)
hAxes.FontSize = fff2;

% Stats
[h,ptotal,ci,stats] = ttest2(ak(1:3), ak(4:6));
cohens_d_md = (mean(ak(1:3)) - mean(ak(4:6))) / stats.sd;
disp(['MD Cohen''s d: ', num2str(cohens_d_md)]);

[h,pvar,ci,stats] = ttest2(vak(1:3), vak(4:6));
cohens_d_rd = (mean(vak(1:3)) - mean(vak(4:6))) / stats.sd;
disp(['RD Cohen''s d: ', num2str(cohens_d_rd)]);

[h,pmicro,ci,stats] = ttest2(uak(1:3), uak(4:6));
cohens_d_ad = (mean(uak(1:3)) - mean(uak(4:6))) / stats.sd;
disp(['AD Cohen''s d: ', num2str(cohens_d_ad)]);
disp(ptotal)
disp(pvar)
disp(pmicro)

adj_p = [ptotal, pvar, pmicro]
%raw_pvals = [ptotal, pvar, pmicro];
%[h_fdr, crit_p, adj_ci_cvrg, adj_p] = fdr_bh(raw_pvals, 0.05, 'pdep', 'no');

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

print(fig, 'figureS3.pdf', '-dpdf', '-painters');
%print(fig, 'figure5.eps', '-depsc', '-painters');
fprintf('Figures saved successfully!\n');

