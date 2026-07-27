clear all
close all
clc

fig = figure('color', [1 1 1], 'Units', 'centimeters', ...
    'Position', [-30 -5 21 25]);
set(fig, 'PaperPositionMode', 'auto');
set(fig, 'PaperOrientation', 'portrait');
set(groot, 'DefaultTextFontName', 'Arial');

fs = filesep;


% Which data
whichd = 1;
dinfo = fun_data_dirs(whichd);

% load data
di = 1;
file_name = dinfo(di).savename;
load([dinfo(1).savename(1:end-6), 'cti_d_c'])
load([dinfo(1).savename(1:end-6), 'MASK_wm_st.mat'])
points_st = points;
load([dinfo(1).savename(1:end-6), 'MASK_wm_ct.mat'])
points_ct = points;

neg = 0;
ssst = 1.5;
sss = 1.5;
fff = 10;

sel_from_mask = squeeze(sum(sum(mask_wm_st, 1), 2));

SIZ = size(FA);

vec = 1:SIZ(3);
sel = vec(sel_from_mask>0);

% Calculating OP and uFA

KANISO = max(KANISO,0); 

numerator = 15.0 * KANISO;
denominator = 10.0 * KANISO + 12.0;

ufa2 = numerator ./ denominator;
ufa2 = max(min(ufa2, 1.0), 0.0);
uFA = sqrt(ufa2);

OP = order_parameter_from_ufa(uFA, FA);
OP(FA>uFA) = 0;
OP(FA<0.2) = 0;

limst = [0, ssst];
limso = [0, sss];

mycmap = turbo(1000);

%% PANEL A
for si=13
    
    ax = subplot_tight(7, 6, 1);
    imagesc(squeeze(MD(:, :, si)), [0 2]); axis image; axis off;
    hold on
    %plot(points_st(si).xi, points_st(si).yi, 'red', 'LineWidth', 2)
    %plot(points_ct(si).xi, points_ct(si).yi, 'red', 'LineWidth', 2)
    title('A1) $$\overline{D} ({\mu}m^2/ms)$$', 'Interpreter','latex', 'fontsize', fff)
    colormap_set
    cb.Ticks = [0 2];
    
    ax = subplot_tight(7, 6, 2);
    imagesc(squeeze(RD(:, :, si)), [0 2]); axis image; axis off;
    hold on
    title({'\textbf{A) CTI maps}','A2) $$D^\bot({\mu}m^2/ms)$$'}, 'Interpreter','latex', 'fontsize', fff)
    colormap_set
    cb.Ticks = [0 2];
    
    ax = subplot_tight(7, 6, 3);
    imagesc(squeeze(AD(:, :, si)), [0 2]); axis image; axis off;
    hold on
    title('A3) $$D^\parallel  ({\mu}m^2/ms)$$', 'Interpreter','latex', 'fontsize', fff)
    colormap_set
    cb.Ticks = [0 2];
    
    ax = subplot_tight(7, 6, 7);
    imagesc(squeeze(MKT(:, :, si)), limso); axis image; axis off;
    hold on
    title('A4) $$\overline{K}_t$$', 'Interpreter','latex', 'fontsize', fff)
    colormap_set
    cb.Ticks = [0 1.5];
    
    ax = subplot_tight(7, 6, 8);
    imagesc(squeeze(RKT(:, :, si)), limso); axis image; axis off;
    hold on
    title('A5) $$K_t^\bot$$', 'Interpreter','latex', 'fontsize', fff)
    colormap_set
    cb.Ticks = [0 1.5];
    
    ax = subplot_tight(7, 6, 9);
    imagesc(squeeze(AKT(:, :, si)), limso); axis image; axis off;
    hold on
    title('A6) $$K_t^\parallel$$', 'Interpreter','latex', 'fontsize', fff)
    colormap_set
    cb.Ticks = [0 1.5];
    
    ax = subplot_tight(7, 6, 13);
    imagesc(squeeze(MKTv(:, :, si)), limso); axis image; axis off;
    hold on
    title('A7) $$\overline{K}_v$$', 'Interpreter','latex', 'fontsize', fff)
    colormap_set
    cb.Ticks = [0 1.5];
    
    ax = subplot_tight(7, 6, 14);
    imagesc(squeeze(RKTv(:, :, si)), limso); axis image; axis off;
    hold on
    title('A8) $$K_v^\bot$$', 'Interpreter','latex', 'fontsize', fff)
    colormap_set
    cb.Ticks = [0 1.5];
    
    ax = subplot_tight(7, 6, 15);
    imagesc(squeeze(AKTv(:, :, si)), limso); axis image; axis off;
    hold on
    title('A9) $$K_v^\parallel$$', 'Interpreter','latex', 'fontsize', fff)
    colormap_set
    cb.Ticks = [0 1.5];
    
    ax = subplot_tight(7, 6, 19);
    imagesc(squeeze(MKTi(:, :, si)), limso); axis image; axis off;
    hold on
    title('A10) $$\overline{K}_{\mu}$$', 'Interpreter','latex', 'fontsize', fff)
    colormap_set
    cb.Ticks = [0 1.5];
    
    ax = subplot_tight(7, 6, 20);
    imagesc(squeeze(RKTi(:, :, si)), limso); axis image; axis off;
    hold on
    title('A11) $$K^\bot_{\mu}$$', 'Interpreter','latex', 'fontsize', fff)
    colormap_set
    cb.Ticks = [0 1.5];
    
    ax = subplot_tight(7, 6, 21);
    imagesc(squeeze(AKTi(:, :, si)), limso); axis image; axis off;
    hold on
    title('A12) $$K^\parallel_{\mu}$$', 'Interpreter','latex', 'fontsize', fff)
    colormap_set
    cb.Ticks = [0 1.5];
    
    ax = subplot_tight(7, 6, 25);
    imagesc(squeeze(FA(:, :, si)), [0 1]); axis image; axis off;
    hold on
    title('A13) $$FA$$', 'Interpreter','latex', 'fontsize', fff)
    colormap_set
    cb.Ticks = [0 1];
    
    ax = subplot_tight(7, 6, 26);
    imagesc(squeeze(uFA(:, :, si)), [0 1]); axis image; axis off;
    hold on
    title('A14) $${\mu}FA$$', 'Interpreter','latex', 'fontsize', fff)
    colormap_set
    cb.Ticks = [0 1];
    
    ax = subplot_tight(7, 6, 27);
    imagesc(squeeze(OP(:, :, si)), [0 1]); axis image; axis off;
    hold on
    title('A15) $$OP (WM)$$', 'Interpreter','latex', 'fontsize', fff)
    colormap_set
    cb.Ticks = [0 1];
    
    colormap(mycmap)
    
    disp(si)
    %pause()
end

%% PANEL B

annotation('textbox', [0.55, 0.9, 0.4, 0.1], ...
    'String', '$$\textbf{B) WM ROI Analysis}$$', ...
    'HorizontalAlignment', 'center', ...
    'BackgroundColor', 'none', ...
    'EdgeColor', 'none', 'Interpreter', 'latex');
xsi = 0.2;
ysi = 0.17;

POS1 = [0.52 0.79 xsi ysi];
POS2 = [0.78 0.79 xsi ysi];
POS3 = [0.52 0.56 xsi ysi];
POS4 = [0.78 0.56 xsi ysi];
POS5 = [0.6 0.33 xsi+0.1 ysi];

rnh_006_final_roi_analysis

set(groot,'defaultAxesTickLabelInterpreter','latex');

c = 8;
x = 1/(2*c+5);
xs = (c+1)*x/2;

ct_color = [0.5 0.5 0.5];
st_color = [1 1 1];

bp = 1:4;
sub_mat = [ones(ns, 1); 2*ones(ns, 1); 3*ones(ns, 1); 4*ones(ns, 1)];

% B1) plot mean kurtosis quatities
hAxes = subplot(4, 4, 3);
hAxes.Position = POS1;

ct_all = [mk_wm(:); kaniso_wm(:); kiso_wm(:); uk_wm(:)];
st_all = [mk_st_wm(:); kaniso_st_wm(:); kiso_st_wm(:); uk_st_wm(:)];
ct = [mmk_wm(:); mkaniso_wm(:); mkiso_wm(:); muk_wm(:)];
st = [mmk_st_wm(:); mkaniso_st_wm(:); mkiso_st_wm(:); muk_st_wm(:)];

hAxes.TickLabelInterpreter = 'latex';
b1 = bar(bp-xs, ct, c*x);
b1.FaceColor = ct_color;
hold on
b2 = bar(bp+xs, st, c*x);
b2.FaceColor = st_color;
plot([sub_mat(:)-xs, sub_mat(:)+xs]', [ct_all(:), st_all(:)]', 'black')
plot(sub_mat(:)-xs, ct_all(:), '.black')
plot(sub_mat(:)+xs, st_all(:), '.black')
xlim([0.5, 4.5])
ylim([0, 4])
set(gca,'XTick', [1 2 3 4], 'YTick',[0 1 2 3]);
xticks([1 2 3 4])
xticklabels({'$$\overline{K}_t$$', ...
    '$$\overline{K}_{ani}$$',...
    '$$\overline{K}_{iso}$$',...
    '$$\overline{K}_{\mu}$$'})
title('B1) Mean Quantities', 'Interpreter', 'latex', ...
    'fontsize', fff)
legend('Contralateral', 'Stroke', 'fontsize', fff)
hAxes.FontSize = fff;

% Stats
[h,ptotal,ci,stats] = ttest(mk_wm(:), mk_st_wm(:));
[h,pvar1,ci,stats] = ttest(kaniso_wm(:), kaniso_st_wm(:));
[h,pvar2,ci,stats] = ttest(kiso_wm(:), kiso_st_wm(:));
[h,pmicro,ci,stats] = ttest(uk_wm(:), uk_st_wm(:));

adj_p = [ptotal, pvar1, pvar2, pmicro]*4
raw_pvals = [ptotal, pvar1, pvar2, pmicro];

[h_fdr, crit_p, adj_ci_cvrg, adj_p] = fdr_bh(raw_pvals, 0.05, 'pdep', 'no');


for mi = 1:4
    if mi == 1
        y_max = max([mk_wm(:); mk_st_wm(:)]) + 0.1; % Adjust 0.1 as needed for spacing
    elseif mi == 2
        y_max = max([kaniso_wm(:); kaniso_st_wm(:)]) + 0.1;
    elseif mi == 3
        y_max = max([kiso_wm(:); kiso_st_wm(:)]) + 0.1;
    elseif mi == 4
        y_max = max([uk_wm(:); uk_st_wm(:)]) + 0.1;
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









%% B2) plot radial kurtosis quatities

c = 6;
x = 1/(2*c+5);
xs = (c+1)*x/2;

bp = 1:3;
sub_mat = [ones(ns, 1); 2*ones(ns, 1); 3*ones(ns, 1)];

hAxes = subplot(4, 4, 4);
hAxes.Position = POS2;

ct_all = [rk_wm(:); vrk_wm(:); urk_wm(:)];
st_all = [rk_st_wm(:); vrk_st_wm(:); urk_st_wm(:)];
ct = [mrk_wm(:); mvrk_wm(:); murk_wm(:)];
st = [mrk_st_wm(:); mvrk_st_wm; murk_st_wm];

hAxes.TickLabelInterpreter = 'latex';
b1 = bar(bp-xs, ct, c*x);
b1.FaceColor = ct_color;
hold on
b2 = bar(bp+xs, st, c*x);
b2.FaceColor = st_color;
plot([sub_mat(:)-xs, sub_mat(:)+xs]', [ct_all(:), st_all(:)]', 'black')
plot(sub_mat(:)-xs, ct_all(:), '.black')
plot(sub_mat(:)+xs, st_all(:), '.black')
xlim([0.5, 3.5])
ylim([0, 4])
set(gca,'XTick', [1 2 3], 'YTick',[0 1 2 3]);
xticks([1 2 3])
xticklabels({'$$K^\bot_t$$', ...
    '$$K^\bot_v$$',...
    '$$K^\bot_{\mu}$$'})
title('B2) Radial Quantities', 'Interpreter', 'latex', ...
    'fontsize', fff)
legend('Contralateral', 'Stroke', 'fontsize', fff)
hAxes.FontSize = fff;

% Stats
[h,ptotal,ci,stats] = ttest(rk_wm(:), rk_st_wm(:));
[h,pvar,ci,stats] = ttest(vrk_wm(:), vrk_st_wm(:));
[h,pmicro,ci,stats] = ttest(urk_wm(:), urk_st_wm(:));
disp(ptotal)
disp(pvar)
disp(pmicro)

adj_p = [ptotal, pvar, pmicro]*3

raw_pvals = [ptotal, pvar, pmicro];

[h_fdr, crit_p, adj_ci_cvrg, adj_p] = fdr_bh(raw_pvals, 0.05, 'pdep', 'no');


for mi = 1:3
    if mi == 1
        y_max = max([rk_wm(:); rk_st_wm(:)]) + 0.1; % Adjust 0.1 as needed for spacing
    elseif mi == 2
        y_max = max([vrk_wm(:); vrk_st_wm(:)]) + 0.1;
    elseif mi == 3
        y_max = max([urk_wm(:); urk_st_wm(:)]) + 0.1;
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






%% B4) plot axial kurtosis quatities
hAxes = subplot(4, 4, 7);
hAxes.Position = POS3;

ct_all = [ak_wm(:); vak_wm(:); uak_wm(:)];
st_all = [ak_st_wm(:); vak_st_wm(:); uak_st_wm(:)];
ct = [mak_wm(:); mvak_wm(:); muak_wm(:)];
st = [mak_st_wm(:); mvak_st_wm; muak_st_wm];

hAxes.TickLabelInterpreter = 'latex';
b1 = bar(bp-xs, ct, c*x);
b1.FaceColor = ct_color;
hold on
b2 = bar(bp+xs, st, c*x);
b2.FaceColor = st_color;
plot([sub_mat(:)-xs, sub_mat(:)+xs]', [ct_all(:), st_all(:)]', 'black')
plot(sub_mat(:)-xs, ct_all(:), '.black')
plot(sub_mat(:)+xs, st_all(:), '.black')
xlim([0.5, 3.5])
ylim([0, 4])
set(gca,'XTick', [1 2 3], 'YTick',[0 1 2 3]);
xticks([1 2 3])
xticklabels({'$$K^\parallel_t$$', ...
    '$$K^\parallel_v$$',...
    '$$K^\parallel_{\mu}$$'})
title('B3) Axial Quantities', 'Interpreter', 'latex', ...
    'fontsize', fff)
legend('Contralateral', 'Stroke', 'fontsize', fff)
hAxes.FontSize = fff;

% Stats
[h,ptotal,ci,stats] = ttest(ak_wm(:), ak_st_wm(:));
[h,pvar,ci,stats] = ttest(vak_wm(:), vak_st_wm(:));
[h,pmicro,ci,stats] = ttest(uak_wm(:), uak_st_wm(:));
disp(ptotal)
disp(pvar)
disp(pmicro)

adj_p = [ptotal, pvar, pmicro]*3

raw_pvals = [ptotal, pvar, pmicro];

[h_fdr, crit_p, adj_ci_cvrg, adj_p] = fdr_bh(raw_pvals, 0.05, 'pdep', 'no');


for mi = 1:3
    if mi == 1
        y_max = max([ak_wm(:); ak_st_wm(:)]) + 0.1; % Adjust 0.1 as needed for spacing
    elseif mi == 2
        y_max = max([vak_wm(:); vak_st_wm(:)]) + 0.1;
    elseif mi == 3
        y_max = max([uak_wm(:); uak_st_wm(:)]) + 0.1;
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


%% PANEL B2
hAxes = subplot(4, 4, 8);
hAxes.Position = POS4;

c = 6;
x = 1/(2*c+5);
xs = (c+1)*x/2;

bp = 1:3;
sub_mat = [ones(ns, 1); 2*ones(ns, 1); 3*ones(ns, 1)];

ct_all = [fa_wm(:); ufa_wm(:); op_wm(:)];
st_all = [fa_st_wm(:); ufa_st_wm(:); op_st_wm(:)];
ct = [mfa_wm(:); mufa_wm(:); mop_wm(:)];
st = [mfa_st_wm(:); mufa_st_wm; mop_st_wm];

hAxes.TickLabelInterpreter = 'latex';
b1 = bar(bp-xs, ct, c*x);
b1.FaceColor = ct_color;
hold on
b2 = bar(bp+xs, st, c*x);
b2.FaceColor = st_color;
plot([sub_mat(:)-xs, sub_mat(:)+xs]', [ct_all(:), st_all(:)]', 'black')
plot(sub_mat(:)-xs, ct_all(:), '.black')
plot(sub_mat(:)+xs, st_all(:), '.black')
xlim([0.5, 3.5])
ylim([0, 1.2])
set(gca,'XTick', [1 2 3], 'YTick',[0 0.3 0.6 0.9]);
xticks([1 2 3])
xticklabels({'$$FA$$', ...
    '$${\mu}FA$$',...
    '$$OP$$'})
title('B4) Anisotropy/Dispersion', 'Interpreter', 'latex', ...
    'fontsize', fff)
legend('Contralateral', 'Stroke', 'fontsize', fff)
hAxes.FontSize = fff;

% Stats
[h,ptotal,ci,stats] = ttest(fa_wm(:), fa_st_wm(:));
[h,pvar,ci,stats] = ttest(ufa_wm(:), ufa_st_wm(:));
[h,pmicro,ci,stats] = ttest(op_wm(:), op_st_wm(:));
disp(ptotal)
disp(pvar)
disp(pmicro)

adj_p = [ptotal, pvar, pmicro]*3

raw_pvals = [ptotal, pvar, pmicro];

[h_fdr, crit_p, adj_ci_cvrg, adj_p] = fdr_bh(raw_pvals, 0.05, 'pdep', 'no');


for mi = 1:3
    if mi == 1
        y_max = max([fa_wm(:); fa_st_wm(:)]) + 0.05; % Adjust 0.1 as needed for spacing
    elseif mi == 2
        y_max = max([ufa_wm(:); ufa_st_wm(:)]) + 0.05;
    elseif mi == 3
        y_max = max([op_wm(:); op_st_wm(:)]) + 0.05;
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






%% B5 axial and radial relative differences

hAxes = subplot(4, 4, 11);
hAxes.Position = POS5;

bp = 1:3;
sub_mat = [ones(ns, 1); 2*ones(ns, 1); 3*ones(ns, 1)];


par_all = [delta_ak(:); delta_vak(:); delta_uak(:)];
per_all = [delta_rk(:); delta_vrk(:); delta_urk(:)];
par = [delta_mak(:); delta_mvak(:); delta_muak(:)];
per = [delta_mrk(:); delta_mvrk(:); delta_murk(:)];

hAxes.TickLabelInterpreter = 'latex';
b1 = bar(bp-xs, per, c*x);
b1.FaceColor = [45 11 117]/255;
hold on
b2 = bar(bp+xs, par, c*x);
b2.FaceColor = [68 154 217]/255;
plot([sub_mat(:)-xs, sub_mat(:)+xs]', [per_all(:), par_all(:)]', 'black')
plot((sub_mat(:)-xs)', (per_all(:))', '.black')
plot((sub_mat(:)+xs)', (par_all(:))', '.black')
xlim([0.5, 3.5])
ylim([0, 3])
xt=sort([bp-xs, bp+xs]);
set(gca,'XTick', xt, 'YTick', [0 1 2 3]);
xticks(xt)
xticklabels({'$$K_t^\bot$$', '$$K_t^\parallel$$',...
    '$$K_v^\bot$$', '$$K_v^\parallel$$',...
    '$$K_{\mu}^\bot$$', '$$K_{\mu}^\parallel$$'})
legend('Radial', 'Axial', 'fontsize', fff, 'location', 'northwest')
title('B5) Metrics Difference', 'Interpreter','latex', 'fontsize', fff)
hAxes.FontSize = fff;

[h,ptotal,ci,stats] = ttest(delta_ak(:), delta_rk(:));
[h,pvar,ci,stats] = ttest(delta_vak(:), delta_vrk(:));
[h,pmicro,ci,stats] = ttest(delta_uak(:), delta_urk(:));
disp(ptotal)
disp(pvar)
disp(pmicro)

adj_p = [ptotal, pvar, pmicro]*3

raw_pvals = [ptotal, pvar, pmicro];

[h_fdr, crit_p, adj_ci_cvrg, adj_p] = fdr_bh(raw_pvals, 0.05, 'pdep', 'no');


for mi = 1:3
    if mi == 1
        y_max = max([delta_ak(:); delta_rk(:)]) + 0.1; % Adjust 0.1 as needed for spacing
    elseif mi == 2
        y_max = max([delta_vak(:); delta_vrk(:)]) + 0.1;
    elseif mi == 3
        y_max = max([delta_uak(:); delta_urk(:)]) + 0.1;
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










%% PANEL C
pppp = 0.25;

annotation('textbox', [0.06, 0.235, 0.4, 0.06], ...
    'String', '$$\textbf{C) Histology - Corpus Callosum}$$', ...
    'HorizontalAlignment', 'center', ...
    'BackgroundColor', 'none', ...
    'EdgeColor', 'none', 'Interpreter', 'latex');

annotation('textbox', [0.52, 0.235, 0.4, 0.06], ...
    'String', '$$\textbf{D) Histology - Barrel Cortex}$$', ...
    'HorizontalAlignment', 'center', ...
    'BackgroundColor', 'none', ...
    'EdgeColor', 'none', 'Interpreter', 'latex');

mainpath = 'C:\Users\rafae\Data\Beading_histo_analysis\';
lw = 2;               % Line width for the arrows
arrows = 0.5;         % Arrow head size
scale_px = 141.5;     % 5 microns = 141.5 pixels
label = '$$5 \mu m$$';
fff = 10;             % Font size for title
margin_x = 0.05;     % 5% margin from the right edge for x-coordinate
margin_y = 0.05;     % 5% margin from the bottom edge for y-coordinate
text_offset_x = -0.1;% Offset for text from the scale bar start (relative to x-range)
text_offset_y = -0.05;% Offset for text from the scale bar y-position (relative to y-range)


%% === First Image: Contralateral Corpus Callosum ===
name = [mainpath, 'Dataset6\ground_truth\images\267.tif.tif'];
t = Tiff(name, 'r');
Image1 = read(t);
ax1 = subplot_tight(3, 4, 9);
POS = ax1.Position;
POS(end) = pppp;
ax1.Position = POS;
imagesc(Image1, [0, 3500])
xlim([500, 3000])
ylim([500, 3000])
axis square
axis off
colormap(ax1, parula)
hold on

% Add scale bar
addScaleBar(scale_px, label, lw, arrows, margin_x, margin_y, text_offset_x+0.05, text_offset_y)

title('C1) Contralateral', 'Interpreter','latex', 'fontsize', fff)

%% === Second Image: Stroke ===
name = [mainpath, 'Dataset5\ground_truth\images\img-221.tif.tif'];
t = Tiff(name, 'r');
Image2 = read(t);

ax1 = subplot_tight(3, 4, 10);
POS = ax1.Position;
POS(end) = pppp;
ax1.Position = POS;
imagesc(Image2, [0, 3500])
xlim([1000, 3500])
ylim([0, 2500])
axis square
axis off
colormap(ax1, parula)
hold on

% Add scale bar
addScaleBar(scale_px, label, lw, arrows, margin_x, margin_y, text_offset_x+0.05, text_offset_y)

title('C2) Stroke', 'Interpreter','latex', 'fontsize', fff)

%% === First Image: Contralateral Cortex ===
for ni = 36 %[11:16, 24:25, 32, 33, 35, 36]
    %for ni = [11:13, 16, 18, 20:26, 29, 31, 32, 34, 35, 37]
    name = [mainpath, 'Dataset7\ground_truth\images\img-',num2str(ni),'.tif.tif'];
    t = Tiff(name, 'r');
    Image1 = read(t);
    
    ax1 = subplot_tight(3, 4, 11);
    POS = ax1.Position;
    POS(end) = pppp;
    ax1.Position = POS;
    imagesc(Image1, [0, 18000])
    xlim([500, 2000])
    ylim([1750, 3250])
    axis square
    axis off
    colormap(ax1, parula)
    hold on
end

% Add scale bar
addScaleBar(scale_px, label, lw, arrows, margin_x, margin_y, text_offset_x, text_offset_y)

title('D1) Contralateral', 'Interpreter', 'latex', 'fontsize', fff)

%%
for ni = 88 %[66, 68, 73, 74, 76:78, 80, 82, 83, 85, 87:90, 92];
    name = [mainpath, 'Dataset1\ground_truth\images\',num2str(ni),'.tif'];
    t = Tiff(name, 'r');
    Image2 = read(t);
    
    ax1 = subplot_tight(3, 4, 12);
    POS = ax1.Position;
    POS(end) = pppp;
    ax1.Position = POS;
    imagesc(Image2, [0, 18000])
    xlim([500, 2000])
    ylim([1000, 2500])
    axis square
    axis off
    colormap(ax1, parula)
    hold on
end

% Add scale bar
addScaleBar(scale_px, label, lw, arrows, margin_x, margin_y, text_offset_x, text_offset_y)

title('D2) Stroke', 'Interpreter', 'latex', 'fontsize', fff)

print(fig, 'figure4.pdf', '-dpdf', '-painters');
%print(fig, 'figure4_final.eps', '-depsc', '-painters');
fprintf('Figures saved successfully!\n');




%% === Helper function: Add a 5 µm scale bar ===
function addScaleBar(length_px, label_text, lw, arrowsize, margin_x, margin_y, text_offset_x, text_offset_y)
% Get current axes limits
x_limits = xlim;
y_limits = ylim;

% Calculate position for the scale bar
% x_start: Position from the right edge, adjusted by margin_x
x_start = x_limits(2) - length_px - (x_limits(2) - x_limits(1)) * margin_x;
% y_pos: Position from the bottom edge, adjusted by margin_y
y_pos = y_limits(2) - (y_limits(2) - y_limits(1)) * margin_y;

% Define the scale bar coordinates
x = [x_start, x_start + length_px];
y = [y_pos, y_pos];

% Calculate text position
text_x = x(1) - (x_limits(2) - x_limits(1)) * text_offset_x;
text_y = y(1) + (y_limits(2) - y_limits(1)) * text_offset_y; % Adjust this for proper vertical alignment

% Add the text label
text(text_x, text_y, label_text, 'BackGroundColor', [1 1 1], 'Interpreter','latex', 'fontsize', 10, 'color', [0 0 0], 'HorizontalAlignment', 'right')

% Plot the scale bar (two quivers for the arrows at each end)
quiver(x(1), y(1), x(2)-x(1), 0, 0, 'color', [0 0 0], 'LineWidth', lw, 'MaxHeadSize', arrowsize)
quiver(x(2), y(2), x(1)-x(2), 0, 0, 'color', [0 0 0], 'LineWidth', lw, 'MaxHeadSize', arrowsize)

end