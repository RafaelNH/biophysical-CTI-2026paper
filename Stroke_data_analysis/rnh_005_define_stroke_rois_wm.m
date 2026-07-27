close all
clear all
clc

addpath(genpath('..\DWIu_v1.11\masktools\'));

fs = filesep;

for whichd = 10
dinfo = fun_data_dirs(whichd);
di = 1;
load([dinfo(di).savename(1:end-6), 'dMASKS.mat'], 'mask_final')
normalize = true;

file_name = [dinfo(di).savename(1:end-6)];

load([file_name, 'dcti_d_c'])

% select diffusion-weighted signals higher b-values
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

figure('color', [1 1 1])
subplot(2, 1, 1)
imagesc(Data_pa_sel(:, :))
subplot(2, 1, 2)
imagesc(Data_pa_sel(:, :)>15)

Data_pa_sel = cat(4, Data_pa_sel, Data_pa_sel);
if whichd == 1
    indexs = [9, 16]; %7
elseif whichd == 2
    indexs = [9, 16]; % 7
elseif whichd == 3
    indexs = [8, 15]; %7
elseif whichd == 4
    indexs = [8, 15]; %7
elseif whichd == 5
    indexs = [9, 15]; % 6
elseif whichd == 6
    indexs = [7, 14]; %7
elseif whichd == 7
    indexs = [8, 14]; %6
elseif whichd == 8
    indexs = [7, 13]; %6
elseif whichd == 9
    indexs = [10, 17]; %7
elseif whichd == 10
    indexs = [9, 16]; %7
end


%% Define initial wm stroke position based on high b-value signals
[mask_wm_st, points] = manual_mask_wpoints(Data_pa_sel, indexs);

points_st = points;

save([dinfo(1).savename(1:end-6), 'MASK_st_lesion.mat'], 'mask_wm_st', 'points')


%% Adjust wm stroke position based on FA maps
SIZ = size(Data_pa_sel);
nd = length(SIZ); % Number of dimensions
ii = 1;
ff = SIZ(3);

fig = figure('color', [1 1 1]);


for sl = indexs(1):indexs(2)
    
    figure(10)
    datas = squeeze(FA(:, :, sl));
    imagesc(datas)
    %colormap('gray')
    axis image
    axis off
    hold on
    pbaspect([1 1 1])
    xi =  points_st(sl).xi;
    yi =  points_st(sl).yi;
    plot(xi, yi, 'red')
    
    [roisi, xi, yi] = roipoly;
    mask_wm_st(:, :, sl) = roisi==1;
    
    plot(xi, yi)
    hold off
    
    points(sl).xi = xi;
    points(sl).yi = yi;
end

save([dinfo(1).savename(1:end-6), 'MASK_wm_st.mat'], 'mask_wm_st', 'points')
points_st = points;

%% Define contralateral area
mask_wm_ct = zeros(SIZ(1:3));

for sl = indexs(1):indexs(2)
    
    figure(10)
    datas = squeeze(FA(:, :, sl));
    imagesc(datas)
    %colormap('gray')
    axis image
    axis off
    hold on
    pbaspect([1 1 1])
    xi =  points_st(sl).xi;
    yi =  points_st(sl).yi;
    plot(xi, yi, 'red')
    
    [roisi, xi, yi] = roipoly;
    mask_wm_ct(:, :, sl) = roisi==1;
    
    plot(xi, yi)
    hold off
    
    points(sl).xi = xi;
    points(sl).yi = yi;
end

save([dinfo(1).savename(1:end-6), 'MASK_wm_ct.mat'], 'mask_wm_ct', 'points')

disp(whichd)


end


