close all
clear all
clc

fs = filesep;
addpath(['..', fs, 'DWIu_v1.11', fs, 'efficient_subpixel_registration'])
addpath(genpath('..\DWIu_v1.11\masktools\'));

denoised = true;
complex = true;
if denoised
    den = '_den';
else
    den = '';
end
if complex
    c = '_complex';
else
    c = '';
end

% Which data
vvv = 100; % numero de direcoes para visualizacao
for whichd = 1
    dinfo = fun_data_dirs(whichd);
    
    produce_mask = true;
    
    
    % sort experiments fromm lower b-values for upper b-values
    maxb = zeros(1, length(dinfo));
    for di = 1:length(dinfo)
        load([dinfo(di).savename, c, den])
        maxb(di) = max(gtab.ub1 + gtab.ub2);
    end
    
    [maxbs, inds] = sort(maxb);
    dinfo = dinfo(inds);
    dinfo(:).savename
    
    tmix = gtab.tmix;
    for di = 1:length(dinfo)
        
        file_name = dinfo(di).savename;
        if di == 1
            load([file_name, c, den])
            disp([file_name, c, den])
            disp(max(gtab.ub1 + gtab.ub2))
            if produce_mask
                [mask_final, points] = manual_mask_wpoints(data);
                if dinfo(di).folder>9
                    save([dinfo(di).savename(1:end-6), 'MASKS.mat'], 'mask_final', 'points')
                else
                    save([dinfo(di).savename(1:end-5), 'MASKS.mat'], 'mask_final', 'points')
                end
            end
            data = align_increasing_b(data, gtab, vvv);
            save([file_name, c, den, '_align'], 'data', 'gtab')
            
        else
            vvv = 100;
            data_ref = data;
            load([file_name, c, den])
            disp([file_name, c, den])
            disp(max(gtab.ub1 + gtab.ub2))
            data = align_to_another(data, data_ref, vvv);
            save([file_name, c, den, '_align'], 'data', 'gtab')
        end
        
    end
end

function rdata = align_increasing_b(data, gtab, vvv)

% Open first dataset
rdata = data;
ivols = 1:size(data, 4);

% parameters to debug
image_velocity = 0.05;

% select lower b-values
ib0s = ivols(gtab.shells==1);
S0s = data(:, :, :, ib0s);
nb0s = size(S0s, 4);
ref = squeeze(S0s(:, :, :, 1));


vieww = 1;
fft2ref = fft2(ref(:, :));
for b0i = 1:nb0s
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
    
    rdata(:, :, :, ib0s(b0i)) = reshape(rimg, size(ref));
end

meanS0s = mean(rdata(:, :, :, ib0s), 4);

refS0 = meanS0s(:, :);
sbi_all = [5, 2, 3, 4];
size_ori = size(meanS0s);

sbii = 1;
for sbi = sbi_all
    
    disp([mean(gtab.bval1(gtab.shells == sbi)),...
        mean(gtab.bval2(gtab.shells == sbi))])
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



