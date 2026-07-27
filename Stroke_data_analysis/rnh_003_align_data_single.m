close all
clear all
clc

fs = filesep;
addpath(['..', fs, 'DWIu_v1.11', fs, 'efficient_subpixel_registration'])
addpath(genpath('..\DWIu_v1.11\masktools\'));

% Which data
vvv = 100;
for whichd = 10
    dinfo = fun_data_dirs(whichd);
    
    produce_mask = true;
    mask_manual = true;
    
    % sort experiments fromm lower b-values for upper b-values
    maxb = zeros(1, length(dinfo));
    for di = 1:length(dinfo)
        load([dinfo(di).savename,'_den'])
        maxb(di) = max(gtab.ub1 + gtab.ub2);
    end
    
    [maxbs, inds] = sort(maxb);
    dinfo = dinfo(inds);
    dinfo(:).savename
    
    tmix = gtab.tmix;
    for di = 1:length(dinfo)
        
        file_name = dinfo(di).savename;
        if di == 1
            load([file_name, '_den'])
            disp([file_name, '_den'])
            disp(max(gtab.ub1 + gtab.ub2))
            if produce_mask
                if mask_manual
                    [mask_final, points] = manual_mask_wpoints(data);
                    save([dinfo(di).savename(1:end-6), 'MASKS.mat'], 'mask_final', 'points')
                else
                    mask_final = simple_thr_mask(data, 0.35);
                    save([dinfo(di).savename(1:end-6), 'MASKS.mat'], 'mask_final')
                    data_sel = data(:, :, :, 1);
                    figure(1);
                    subplot(411); imagesc(data_sel(:, :)); axis image; axis off; colormap(gray(256));
                    subplot(412); imagesc(mask_final(:, :)); axis image; axis off; colormap(gray(256));
                end
            end
            data = align_increasing_b(data, gtab, vvv);
            save([file_name, '_align'], 'data', 'gtab')
            
        else
            data_ref = data;
            load([file_name, '_den'])
            disp([file_name, '_den'])
            disp(max(gtab.ub1 + gtab.ub2))
            data = align_to_another(data, data_ref, vvv);
            save([file_name, '_align'], 'data', 'gtab')
        end
        
    end
end

function rdata = align_increasing_b(data, gtab, vvv)

% Open first dataset
rdata = data;
ivols = 1:size(data, 4);
%N = size(data, 3);

% parameters to debug
image_velocity = 0.05;

% first lower b-value
[bm, ind] = min(gtab.bval);
refb0 = squeeze(data(:, :, :, ind));

% select lower b-values (should be b-value=0)
ib0s = ivols(gtab.bval == bm);
S0s = data(:, :, :, ib0s);
nb0s = size(S0s, 4);

% align lower b-values (should be b-value=0)
%[optimizer, metric] = imregconfig('monomodal');
%for si = 1:N
vieww = 1;
%ref = refb0(:, :, si);
ref = refb0;
fft2ref = fft2(ref(:, :));
for b0i = 2:nb0s
    %cimg = squeeze(S0s(:, :, si, b0i));
    cimg = squeeze(S0s(:, :, :, b0i));
    cimg = cimg(:, :);
    [output, nimg] = dftregistration(fft2ref, fft2(cimg), 10);
    rimg = abs(ifft2(nimg));
    %rimg = imregister(cimg, ref, 'translation', optimizer, metric);
    
    if vieww < vvv
        figure(1);
        subplot(411); imagesc(cimg); axis image; axis off; title(['original ', num2str(b0i)]); colormap(gray(256));
        subplot(412); imagesc(rimg); axis image; axis off; title('corrected'); colormap(gray(256));
        subplot(413); imagesc(cimg-ref(:, :), [-0.4 0.4]); axis image; axis off; title('ori-ref'); colormap(gray(256));
        subplot(414); imagesc(cimg-rimg); axis image; axis off; title('ori-alig'); colormap(gray(256));
        %colorbar
        drawnow; pause(image_velocity)
        vieww = vieww + 1;
    end
    
    rdata(:, :, :, ib0s(b0i)) = reshape(rimg, size(refb0));
end
%end

meanS0s = mean(rdata(:, :, :, ib0s), 4);

refS0 = meanS0s(:, :);
[sb, sbi_all] = sort(gtab.ub1+gtab.ub2);
size_ori = size(meanS0s);

sbii = 1;
for sbi = sbi_all
    
    idwi = ivols(gtab.shells == sbi);
    dwi_sel = data(:, :, :, idwi);
    ndwi = size(dwi_sel, 4);
    
    vieww = 1;
    if sbii == 1
        ref = refS0;
        fft2ref = fft2(refS0);
    end
    for bi = 1:ndwi
        if sbii > 1
            ref = squeeze(dwi_ref(:, :, :, bi));
            fft2ref = fft2(ref(:, :));
        end
        cimg = squeeze(dwi_sel(:, :, :, bi));
        cimg = cimg(:, :);
        [output, nimg] = dftregistration(fft2ref, fft2(cimg), 10);
        rimg = abs(ifft2(nimg));
        %rimg = imregister(cimg, ref, 'translation', optimizer, metric);
        
        if vieww < vvv
            figure(1);
            subplot(411); imagesc(cimg); axis image; axis off; title(['original ', num2str(bi)]); colormap(gray(256));
            subplot(412); imagesc(rimg); axis image; axis off; title('corrected'); colormap(gray(256));
            subplot(413); imagesc(rimg-ref(:, :)); axis image; axis off; title('alig-ref'); colormap(gray(256));
            subplot(414); imagesc(cimg-rimg); axis image; axis off; title('ori-alig'); colormap(gray(256));
            %colorbar
            drawnow; pause(image_velocity)
            vieww = vieww + 1;
        end
        
        rdata(:, :, :, idwi(bi)) = reshape(rimg, size_ori);
        
    end
    sbii = sbii + 1;
    dwi_ref = rdata(:, :, :, idwi);
end
end



function rdata = align_to_another(data, data_ref, vvv)

% Open first dataset
rdata = data;
nvols = size(data, 4);
size_ori = size(data);

% parameters to debug
image_velocity = 0.1;
vieww = 1;
for vi = 1:nvols
    ref = squeeze(data_ref(:, :, :, vi));
    fft2ref = fft2(ref(:, :));
    
    cimg = squeeze(data(:, :, :, vi));
    cimg = cimg(:, :);
    [output, nimg] = dftregistration(fft2ref, fft2(cimg), 10);
    rimg = abs(ifft2(nimg));
    if vieww < vvv
        figure(1);
        subplot(411); imagesc(cimg); axis image; axis off; title(['original ', num2str(vi)]); colormap(gray(256));
        subplot(412); imagesc(rimg); axis image; axis off; title('corrected'); colormap(gray(256));
        subplot(413); imagesc(rimg-ref(:, :)); axis image; axis off; title('alig-ref'); colormap(gray(256));
        subplot(414); imagesc(cimg-rimg); axis image; axis off; title('ori-alig'); colormap(gray(256));
        %colorbar
        drawnow; pause(image_velocity)
        vieww = vieww + 1;
    end
    
    rdata(:, :, :, vi) = reshape(rimg, size_ori(1:3));
    
end
end



