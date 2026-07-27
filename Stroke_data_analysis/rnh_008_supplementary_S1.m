close all
clear all
clc

addpath ../DWIu_toolbox_v1.11/DKI_toolbox/
addpath ../DWIu_toolbox_v1.11/DDE_toolbox/

fff = 10;
fff2 = 10;
fig = figure('color', [1 1 1], 'Units', 'centimeters', ...
    'Position', [1 1 21 6]);

set(fig, 'PaperPositionMode', 'auto');
set(fig, 'PaperOrientation', 'portrait');
set(groot, 'DefaultTextFontName', 'Arial');

fs = filesep;


% Which data
whichd = 1;
dinfo = fun_data_dirs(whichd);

% load data
for di = 1:length(dinfo)
    file_name = [dinfo(di).savename];
    pt = [file_name, '_align'];
    
    load(pt)
    
    sized = size(data);
    nb = length([0, gtab.ub1]);
    sized(4) = nb;
    
    if di == 1
        b1_all = [0, gtab.ub1];
        
        Data_pa_nonorm = zeros(sized);
        for bi = 0:(nb-1)
            Data_pa_nonorm(:, :, :, bi+1) = mean(data(:, :, :, gtab.shells==bi), 4);
        end
        
    else
        
        b1_new = [0, gtab.ub1];
        Data_pa_nonorm_new = zeros(sized);
        for bi = 0:(nb-1)
            Data_pa_nonorm_new(:, :, :, bi+1) = mean(data(:, :, :, gtab.shells==bi), 4);
        end
        
        b1_all = [b1_all, b1_new];
        Data_pa_nonorm = cat(4, Data_pa_nonorm, Data_pa_nonorm_new);
        
    end
    
end

[mbval, maxi] = max(b1_all);

Data_pa_sel = Data_pa_nonorm(:, :, :, maxi);

di = 1;
file_name = dinfo(di).savename;
load([dinfo(1).savename(1:end-6), 'cti_d_c'])
load([dinfo(1).savename(1:end-6), 'MASK_st_lesion.mat'], 'mask_wm_st', 'points')
points_lesion = points;
load([dinfo(1).savename(1:end-6), 'MASK_wm_st.mat'])
points_st = points;
load([dinfo(1).savename(1:end-6), 'MASK_wm_ct.mat'])
points_ct = points;




%

neg = 0;
ssst = 1.5;
sss = 1.5;
fff = 10;

sel_from_mask = squeeze(sum(sum(mask_wm_st, 1), 2));

SIZ = size(FA);

vec = 1:SIZ(3);
sel = vec(sel_from_mask>0);


limst = [0, ssst];
limso = [0, sss];

mycmap = turbo(1000);

for si=13
    ax = subplot_tight(1, 5, 1);
    imagesc(data(:, :, si, 1)); axis image; axis off;
    title({'A) $$S(0, 0)$$'}, 'Interpreter','latex', 'fontsize', fff)
    colormap(ax, gray)
    
    ax = subplot_tight(1, 5, 2);
    imagesc(Data_pa_sel(:, :, si)); axis image; axis off;
    title({'B) $$\overline{S}(3, 0)$$'}, 'Interpreter','latex', 'fontsize', fff)
    colormap(ax, gray)
    hold on
    plot(points_lesion(si).xi, points_lesion(si).yi, 'blue', 'LineWidth', 1)
    
    ax = subplot_tight(1, 5, 3);
    imagesc(squeeze(MD(:, :, si)), [0 1.5]); axis image; axis off;
    hold on
    plot(points_lesion(si).xi, points_lesion(si).yi, 'blue', 'LineWidth', 1)
    title({'C) $$\overline{D}({\mu}m^2/ms)$$'}, 'Interpreter','latex', 'fontsize', fff)
    colormap(ax, gray)
    pos = get(ax, 'position');
    cb = colorbar('southoutside');
    cb.Position =  [pos(1) pos(2)+0.15 pos(3) 0.05];
    cb.FontSize = fff;
    cb.Ticks = [0 1.5];
    
    ax = subplot_tight(1, 5, 4);
    imagesc(squeeze(FA(:, :, si)), [0 1]); axis image; axis off;
    hold on
    plot(points_lesion(si).xi, points_lesion(si).yi, 'blue', 'LineWidth', 1)
    plot(points_st(si).xi, points_st(si).yi, 'red', 'LineWidth', 1)
    plot(points_ct(si).xi, points_ct(si).yi, 'green', 'LineWidth', 1)
    title('D) $$FA$$', 'Interpreter','latex', 'fontsize', fff)
    colormap(ax, gray)
    pos = get(ax, 'position');
    cb = colorbar('southoutside');
    cb.Position =  [pos(1) pos(2)+0.15 pos(3) 0.05];
    cb.FontSize = fff;
    cb.Ticks = [0 1];
    lg = legend('Stroke Lesion','WM (Stroke)', 'WM (Contralateral)', ...
        'location','eastoutside');
    pos = lg.Position;
    lg.Position = [pos(1)+0.15 pos(2)+0.05 pos(3) pos(4)];
   
end

print(fig, 'figureS1.pdf', '-dpdf', '-painters');
%print(fig, 'figureS1.eps', '-depsc', '-painters');
fprintf('Figures saved successfully!\n');
