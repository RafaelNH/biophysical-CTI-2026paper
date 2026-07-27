close all
clc

% this code was run before some changes in the organization of the data
% some directories may have to be addapted

addpath('..\DWIu_v1.9\DDE_toolbox')
addpath(genpath('..\DWIu_v1.9\pvtools/'));
addpath(genpath('..\DWIu_v1.9\masktools/'));

% Which data
whichd = 6;
[folder_name, folders, EXPS, refEXPS, outname] = fun_data_dirs(whichd);

%load([folder_name, filesep, 'MASKS.mat'], 'mask_final')
ali = 'align_';
load(['cti_', ali, outname])

siz = size(KTOTAL_pa);
KTOTAL_pa_hack = zeros([siz, 2]);
KTOTAL_pa_hack(:, :, :, 1) = KTOTAL_pa;

%% optional - define mask and regions of interest
[tumour_mask, points] = manual_mask_wpoints(KTOTAL_pa_hack);
save([folder_name, filesep, 'tumour_mask.mat'], 'tumour_mask', 'points')

%[gm_mask, points] = manual_mask_wpoints(KTOTAL_pa_hack);
%save([folder_name, filesep, 'gm_mask.mat'], 'gm_mask', 'points')


function [mask, mask_points] = manual_mask_wpoints(image, indexs)
% Mask Manual delination
%
%     Parameters
%     ----------
%     image : array (..., N)
%         Array containing the diffusion-weighted signals, note that last
%         dimension corresponds to different diffusion-weigthed experiment
%         which may have been acquired for different gradient direction
%         or b-value
%     indexs : array [ii, ff], optional
%         Array containing the initial and final slice indexes
%         Default: All slices
%
%     Returns
%     -------
%     mask : array (..., )
%         Array containing true values for voxels to be processed
%

SIZ = size(image);
nd = length(SIZ); % Number of dimensions

fig = figure('color', [1 1 1]);
if nd == 4
    if ~exist('indexs', 'var')
        ii = 1;
        ff = SIZ(3);
    else
        ii = indexs(1);
        ff = indexs(2);
    end
    
    mask = zeros(SIZ(1:3));
    
    for sl = ii:ff
        
        figure(fig)
        datas = squeeze(image(:, :, sl, 1));
        imagesc(datas, [0 1.5])
        %colormap('gray')
        axis image
        axis off
        hold on
        pbaspect([1 1 1])
        
        [roisi, xi, yi] = roipoly;
        mask(:, :, sl) = roisi==1;
        
        plot(xi, yi)
        hold off
        
        mask_points(sl).xi = xi;
        mask_points(sl).yi = yi;
    end
    
elseif nd == 3
    
    datas = squeeze(image(:, :, 1));
    
    imagesc(datas)
    colormap('gray')
    axis image
    axis off
    hold on
    
    pbaspect([1 1 1])
    
    [roisi, xi, yi] = roipoly;
    
    mask = roisi==1;
    plot(xi, yi)
    
    hold off
    
elseif nd == 2
    
    datas = image;
    imagesc(datas)
    colormap('gray')
    axis image
    axis off
    hold on
    pbaspect([1 1 1])
    
    [roisi, xi, yi] = roipoly;
    
    mask = roisi==1;
    
    plot(xi, yi, 'red')
    
    hold off
    
    mask_points.xi = xi;
    mask_points.yi = yi;
end
end