close all
clear all
clc

fs = filesep;
addpath(genpath(['..', fs, 'DWIu_v1.11', fs, 'pvtools/']));
addpath(['..', fs, 'DWIu_v1.11', fs, 'DDE_toolbox'])
addpath(['..', fs, 'DWIu_v1.11', fs, 'DKI_toolbox'])
addpath(['..', fs, '..', fs, 'diffusion_mri', fs, 'reconst'])
addpath(['..', fs, '..', fs, 'diffusion_mri', fs, 'dirs'])
load('dirs45.mat')

slice = 10;
constraints = true;
Smooth_Gaussian = false;
exlude_b0s = true;
denoised_data = true;

if constraints
    cst = '_c';
else
    cst = '';
end

if Smooth_Gaussian
    gst = '_gs';
else
    gst = '';
end

if denoised_data == true
    den = '_den';
else
    den = '';
end

for whichd = [1:6, 8]
    dinfo = fun_data_dirs(whichd);
    if whichd == 1
        load([dinfo(1).savename(1:end-7), 'MASKS.mat'], 'mask_final')
    elseif whichd == 4
        load([dinfo(1).savename(1:end-5), 'MASKS.mat'], 'mask_final')
    else
        load([dinfo(1).savename(1:end-6), 'MASKS.mat'], 'mask_final')
    end
    normalize = true;
    
    for di = 1:length(dinfo)
        file_name = dinfo(di).savename;
        pt = [file_name, den, '_align'];
        disp(pt)
        
        load(pt)
        if Smooth_Gaussian
            fwhm = 1.25;
            sd = fwhm/(sqrt(8*log(2))); % Convert FWHM to sd
            S = size(data);
            Nvol = S(4);
            size_k = [7 7 7];
            for im = Nvol:-1:1
                for si = S(3):-1:1
                    A = squeeze(data(:,:,si,im));
                    W = imgaussfilt(A, sd);
                    data(:,:,si,im) = W;
                end
            end
        end
        sized = size(data);
        nb = length([0, gtab.ub1]);
        sized(4) = nb;
        
        if di == 1
            bvals1_all = gtab.bval1;
            bvals2_all = gtab.bval2;
            dir1_all = gtab.dir1;
            dir2_all = gtab.dir2;
            b1_all = [0, gtab.ub1];
            b2_all = [0, gtab.ub2];
            cos2ang = [1, gtab.uc2a];
            
            Data_pa_nonorm = zeros(sized);
            for bi = 0:(nb-1)
                Data_pa_nonorm(:, :, :, bi+1) = mean(data(:, :, :, gtab.shells==bi), 4);
            end
            
            mS0 = mean(data(:, :, :, gtab.bval==0), 4);
            for vi = 1:size(data, 4)
                data(:, :, :, vi) = data(:, :, :, vi) ./ mS0;
            end
            Data_all = data;
            
            
            Data_pa = zeros(sized);
            for bi = 0:(nb-1)
                Data_pa(:, :, :, bi+1) = mean(data(:, :, :, gtab.shells==bi), 4);
            end
            
            %             figure
            %             for bi = 1:nb
            %                 subplot(nb, 1, bi)
            %                 data_si = Data_pa(:, :, :, bi);
            %                 imagesc(data_si(:, :), [0 1.2])
            %             end
            
            
            dd = diag(dir1_all'*dir2_all);
            ddd = dd((bvals1_all + bvals2_all)> 0);
            
            %figure, plot(bvals1_all), hold on, plot(bvals2_all, '--'), plot(dd - 2)
            
        else
            
            bvals1_new = gtab.bval1;
            bvals2_new = gtab.bval2;
            dir1_new = gtab.dir1;
            dir2_new = gtab.dir2;
            b1_new = [0, gtab.ub1];
            b2_new = [0, gtab.ub2];
            cos2ang_new = [1, gtab.uc2a];
            
            Data_pa_nonorm_new = zeros(sized);
            for bi = 0:(nb-1)
                Data_pa_nonorm_new(:, :, :, bi+1) = mean(data(:, :, :, gtab.shells==bi), 4);
            end
            
            mS0 = mean(data(:, :, :, gtab.bval==0), 4);
            for vi = 1:size(data, 4)
                data(:, :, :, vi) = data(:, :, :, vi) ./ mS0;
            end
            
            sized = size(data);
            nb = length(b1_new);
            sized(4) = nb;
            Data_pa_new = zeros(sized);
            for bi = 0:(nb-1)
                Data_pa_new(:, :, :, bi+1) = mean(data(:, :, :, gtab.shells==bi), 4);
            end
            
            
            
            %             figure
            %             for bi = 1:nb
            %                 subplot(nb, 1, bi)
            %                 data_si = Data_pa_new(:, :, :, bi);
            %                 imagesc(data_si(:, :), [0 1.2])
            %             end
            
            
            dd = diag(dir1_new'*dir2_new);
            ddd = dd((bvals1_new + bvals2_new)> 0);
            
            %figure, plot(bvals1_new), hold on, plot(bvals2_new, '--'), plot(dd - 2)
            
            
            Data_all = cat(4, Data_all, data);
            bvals1_all = [bvals1_all, bvals1_new];
            bvals2_all = [bvals2_all, bvals2_new];
            dir1_all = [dir1_all, dir1_new];
            dir2_all = [dir2_all, dir2_new];
            b1_all = [b1_all, b1_new];
            b2_all = [b2_all, b2_new];
            cos2ang = [cos2ang, cos2ang_new];
            Data_pa = cat(4, Data_pa, Data_pa_new);
            Data_pa_nonorm = cat(4, Data_pa_nonorm, Data_pa_nonorm_new);
            
        end
        
    end
    
    dd = diag(dir1_all'*dir2_all);
    ddd = dd((bvals1_all + bvals2_all)> 0);
    
    figure, plot(bvals1_all), hold on, plot(bvals2_all, '--'), plot(dd - 2)
    
    if exlude_b0s
        include_d = (bvals1_all + bvals2_all)> 0;
        Data_all = Data_all(:, :, :, include_d==1);
        bvals1_all = bvals1_all(include_d==1);
        bvals2_all = bvals2_all(include_d==1);
        dir1_all = dir1_all(:, include_d==1);
        dir2_all = dir2_all(:, include_d==1);
    end
    if constraints
        [D11,D22,D33,D12,D13,D23,...
            W1111, W2222, W3333, W1112, W1113,...
            W1222, W2223, W1333, W2333, W1122,...
            W1133, W2233, W1123, W1223, W1233,...
            K1111, K2222, K3333, K1112, K1113, ...
            K1222, K2223, K1333, K2333, K1122, ...
            K1133, K2233, K1123, K1322, K1233, ...
            K1212, K1313, K2323, K1223, K1323, K1213,S0, SSE, Ndiff, Nkurt, Ncov]=...
            fun_DDE_CLLS_comp_b2_improved_rh(Data_all, mask_final, ...
            bvals1_all', bvals2_all', dir1_all, dir2_all, V);
    else
        [D11,D22,D33,D12,D13,D23,...
            W1111, W2222, W3333, W1112, W1113,...
            W1222, W2223, W1333, W2333, W1122,...
            W1133, W2233, W1123, W1223, W1233,...
            K1111, K2222, K3333, K1112, K1113, ...
            K1222, K2223, K1333, K2333, K1122, ...
            K1133, K2233, K1123, K1322, K1233, ...
            K1212, K1313, K2323, K1223, K1323, K1213, S0, SSE] = ...
            fun_DDE_ULLS_comp_b2_rh(Data_all, mask_final, ...
            bvals1_all', bvals2_all', dir1_all, dir2_all);
    end
    
    N = size(Data_all);
    
    pars0_all = zeros(N(1), N(2), N(3), 43);
    pars0_all(:, :, :, 1) = S0;
    pars0_all(:, :, :, 2) = D11;
    pars0_all(:, :, :, 3) = D22;
    pars0_all(:, :, :, 4) = D33;
    pars0_all(:, :, :, 5) = D12;
    pars0_all(:, :, :, 6) = D13;
    pars0_all(:, :, :, 7) = D23;
    pars0_all(:, :, :, 8) = W1111;
    pars0_all(:, :, :, 9) = W2222;
    pars0_all(:, :, :, 10) = W3333;
    pars0_all(:, :, :, 11) = W1112;
    pars0_all(:, :, :, 12) = W1113;
    pars0_all(:, :, :, 13) = W1222;
    pars0_all(:, :, :, 14) = W2223;
    pars0_all(:, :, :, 15) = W1333;
    pars0_all(:, :, :, 16) = W2333;
    pars0_all(:, :, :, 17) = W1122;
    pars0_all(:, :, :, 18) = W1133;
    pars0_all(:, :, :, 19) = W2233;
    pars0_all(:, :, :, 20) = W1123;
    pars0_all(:, :, :, 21) = W1223;
    pars0_all(:, :, :, 22) = W1233;
    pars0_all(:, :, :, 23) = K1111;
    pars0_all(:, :, :, 24) = K2222;
    pars0_all(:, :, :, 25) = K3333;
    pars0_all(:, :, :, 26) = K1112;
    pars0_all(:, :, :, 27) = K1113;
    pars0_all(:, :, :, 28) = K1222;
    pars0_all(:, :, :, 29) = K2223;
    pars0_all(:, :, :, 30) = K1333;
    pars0_all(:, :, :, 31) = K2333;
    pars0_all(:, :, :, 32) = K1122;
    pars0_all(:, :, :, 33) = K1133;
    pars0_all(:, :, :, 34) = K2233;
    pars0_all(:, :, :, 35) = K1123;
    pars0_all(:, :, :, 36) = K1322;
    pars0_all(:, :, :, 37) = K1233;
    pars0_all(:, :, :, 38) = K1212;
    pars0_all(:, :, :, 39) = K1313;
    pars0_all(:, :, :, 40) = K2323;
    pars0_all(:, :, :, 41) = K1223;
    pars0_all(:, :, :, 42) = K1323;
    pars0_all(:, :, :, 43) = K1213;
    %pars0_all = real(pars0_all);
    
    
    
    
    [KTOTAL, KANISO, KISO, KINTRA, Conv, MD] = fun_resolve_kurtosis(pars0_all);
    scal = (length(parula):-1:1)/length(parula);
    
    %mycmap = [scal', zeros(length(parula), 2); 0, 0, 0; parula];
    
    mycmap = [0, 0, 0; parula];
    
    figure('color', [1 1 1])
    subplot(4, 1, 1)
    imagesc(squeeze(KTOTAL(:, :, slice)), [0 1.5]); axis image; axis off;
    colormap(mycmap)
    colorbar
    
    subplot(4, 1, 2)
    imagesc(squeeze(KANISO(:, :, slice)), [0 1.5]); axis image; axis off;
    colormap(mycmap)
    colorbar
    
    subplot(4, 1, 3)
    imagesc(squeeze(KISO(:, :, slice)), [0 1.5]); axis image; axis off;
    colormap(mycmap)
    colorbar
    
    subplot(4, 1, 4)
    imagesc(squeeze(KINTRA(:, :, slice)), [0 1.5]); axis image; axis off;
    colormap(mycmap)
    colorbar
    
    if constraints
        
        figure('color', [1 1 1])
        subplot(4, 1, 1)
        imagesc(squeeze(3*45 - Ndiff(:, :, slice)), [0 10]); axis image; axis off;
        colormap(mycmap)
        colorbar
        
        subplot(4, 1, 2)
        imagesc(squeeze(3*45 - Nkurt(:, :, slice)), [0 10]); axis image; axis off;
        colormap(mycmap)
        colorbar
        
        subplot(4, 1, 3)
        imagesc(squeeze(1 - Ncov(:, :, slice)), [0 10]); axis image; axis off;
        colormap(mycmap)
        colorbar
        
    end
    [MD, AD, RD, FA, V1] = dti_metrics(pars0_all(:, :, :, 2:7), mask_final);
    
    MDs = MD(:, :, :);
    ADs = AD(:, :, :);
    RDs = RD(:, :, :);
    FAs = FA(:, :, :);
    
    figure('color', [1 1 1])
    subplot(4, 1, 1)
    imagesc(squeeze(MDs(:, :, slice)), [0 1.5]); axis image; axis off;
    title('MD')
    colorbar
    
    subplot(4, 1, 2)
    imagesc(squeeze(ADs(: ,:, slice)), [0 1.5]); axis image; axis off;
    title('AD')
    colorbar
    
    subplot(4, 1, 3)
    imagesc(squeeze(RDs(:, :, slice)), [0 1.5]); axis image; axis off;
    title('RD')
    colorbar
    
    subplot(4, 1, 4)
    imagesc(squeeze(FAs(:, :, slice)), [0 1]); axis image; axis off;
    title('FA')
    colorbar
    
    colormap(gray)
    
    
    %%
    % Vars = 5/6 * KANISOs .* MDs .* MDs;
    %
    % uFA = sqrt(3/2 * (5/6 * KANISOs)./ ((5/6 * KANISOs) + 1) );
    %
    % figure('color', [1 1 1])
    % subplot(2, 1, 1)
    % imagesc(squeeze(Vars(:, :)), [0 4]); axis image; axis off;
    % title('Vars')
    %
    %
    % subplot(2, 1, 2)
    % imagesc(squeeze(real(uFA(:, :))), [0 1]); axis image; axis off;
    % title('FA')
    % colorbar
    
    
    
    %%
    
    [DT, KT, KTi, KTv] = fun_resolve_kurtosis_tensors(pars0_all);
    
    [MK, AK, RK] = dki_metrics(DT, KT, mask_final, [-1.5, 10]);
    [MKv, AKv, RKv] = dki_metrics(DT, KTv, mask_final, [-1.5, 10]);
    [MKi, AKi, RKi] = dki_metrics(DT, KTi, mask_final, [-1.5, 10]);
    
    [MKT, AKT, RKT] = dki_tensor_metrics(DT, KT, mask_final);
    [MKTv, AKTv, RKTv] = dki_tensor_metrics(DT, KTv, mask_final);
    [MKTi, AKTi, RKTi] = dki_tensor_metrics(DT, KTi, mask_final);
    
    
    sss = 1.5;
    
    scal = (length(parula):-1:1)/length(parula);
    mycmap = [scal', zeros(length(parula), 2); 0, 0, 0; parula];
    
    figure('color', [1 1 1])
    subplot(3, 1, 1)
    imagesc(squeeze(RK(:, :, slice)), [-sss sss]); axis image; axis off;
    title('RK')
    
    subplot(3, 1, 2)
    imagesc(squeeze(AK(:, :, slice)), [-sss sss]); axis image; axis off;
    title('AK')
    
    subplot(3, 1, 3)
    imagesc(squeeze(MK(:, :, slice)), [-sss sss]); axis image; axis off;
    title('MK')
    colormap(mycmap)
    
    
    figure('color', [1 1 1])
    subplot(3, 1, 1)
    imagesc(squeeze(RKT(:, :, slice)), [-sss sss]); axis image; axis off;
    title('RKT')
    
    subplot(3, 1, 2)
    imagesc(squeeze(AKT(:, :, slice)), [-sss sss]); axis image; axis off;
    title('AKT')
    
    subplot(3, 1, 3)
    imagesc(squeeze(MKT(:, :, slice)), [-sss sss]); axis image; axis off;
    title('MKT')
    colormap(mycmap)
    
    
    
    figure('color', [1 1 1])
    subplot(3, 1, 1)
    imagesc(squeeze(RKv(:, :,  slice)), [-sss sss]); axis image; axis off;
    title('RKv')
    
    subplot(3, 1, 2)
    imagesc(squeeze(AKv(:, :, slice)), [-sss sss]); axis image; axis off;
    title('AKv')
    
    subplot(3, 1, 3)
    imagesc(squeeze(MKv(:, :, slice)), [-sss sss]); axis image; axis off;
    title('MKv')
    colormap(mycmap)
    
    
    
    figure('color', [1 1 1])
    subplot(3, 1, 1)
    imagesc(squeeze(RKTv(:, :,  slice)), [-sss sss]); axis image; axis off;
    title('RKTv')
    
    subplot(3, 1, 2)
    imagesc(squeeze(AKTv(:, :, slice)), [-sss sss]); axis image; axis off;
    title('AKTv')
    
    subplot(3, 1, 3)
    imagesc(squeeze(MKTv(:, :, slice)), [-sss sss]); axis image; axis off;
    title('MKTv')
    colormap(mycmap)
    
    
    
    figure('color', [1 1 1])
    subplot(3, 1, 1)
    imagesc(squeeze(RKi(:, :, slice)), [-sss sss]); axis image; axis off;
    title('RKi')
    
    subplot(3, 1, 2)
    imagesc(squeeze(AKi(:, :, slice)), [-sss sss]); axis image; axis off;
    title('AKi')
    
    subplot(3, 1, 3)
    imagesc(squeeze(MKi(:, :, slice)), [-sss sss]); axis image; axis off;
    title('MKi')
    colormap(mycmap)
    
    %for slice = 1:16
    figure('color', [1 1 1])
    subplot(3, 1, 1)
    imagesc(squeeze(RKTi(:, :, slice)), [-sss sss]); axis image; axis off;
    title('RKTi')
    
    subplot(3, 1, 2)
    imagesc(squeeze(AKTi(:, :, slice)), [-sss sss]); axis image; axis off;
    title('AKTi')
    
    subplot(3, 1, 3)
    imagesc(squeeze(MKTi(:, :, slice)), [-sss sss]); axis image; axis off;
    title('MKTi')
    colormap(mycmap)%, end
    
    %% RICE
    %[D2mat, S0mat, S2mat, S4mat, A0mat, A2mat, SSV, COVi] = fun_resolve_invariants(pars0_all, mask_final);
    [D0, D2, vW0, vW2, vW4, uW0, uW2, uW4, A0, A2, At, Cov_ss_n] = fun_resolve_cti_invariants(pars0_all, mask_final);
    
    
    figure('color', [1 1 1])
    subplot(3, 1, 1)
    imagesc(squeeze(D0(:, :, slice)), [0 1.5]); axis image; axis off;
    colorbar
    
    subplot(3, 1, 2)
    imagesc(squeeze(D2(:, :, slice)), [0 1.5]); axis image; axis off;
    colorbar
    
    figure('color', [1 1 1])
    subplot(3, 1, 1)
    imagesc(squeeze(vW0(:, :, slice)), [-sss sss]); axis image; axis off;
    title('vW0')
    colorbar
    
    subplot(3, 1, 2)
    imagesc(squeeze(vW2(:, :, slice)), [-sss/10 sss/10]); axis image; axis off;
    title('vW2')
    colorbar
    
    subplot(3, 1, 3)
    imagesc(squeeze(vW4(:, :, slice)), [-sss/10 sss/10]); axis image; axis off;
    title('vW4')
    colormap(mycmap)
    colorbar
    
    figure('color', [1 1 1])
    subplot(3, 1, 1)
    imagesc(squeeze(uW0(:, :, slice)), [-sss sss]); axis image; axis off;
    title('vW0')
    colorbar
    
    subplot(3, 1, 2)
    imagesc(squeeze(uW2(:, :, slice)), [-sss/10 sss/10]); axis image; axis off;
    title('vW2')
    colorbar
    
    
    subplot(3, 1, 3)
    imagesc(squeeze(uW4(:, :, slice)), [-sss/10 sss/10]); axis image; axis off;
    title('vW4')
    colormap(mycmap)
    colorbar
    
    figure('color', [1 1 1])
    subplot(3, 1, 1)
    imagesc(squeeze(A0(:, :, slice)), [-sss sss]); axis image; axis off;
    title('A0')
    colormap(mycmap)
    colorbar
    
    
    subplot(3, 1, 2)
    imagesc(squeeze(A2(:, :, slice)), [-sss sss]); axis image; axis off;
    colormap(mycmap)
    title('A2')
    colorbar
    
    subplot(3, 1, 3)
    imagesc(squeeze(At(:, :, slice)), [-sss sss]); axis image; axis off;
    colormap(mycmap)
    title('At')
    colorbar
    
    Cov = Cov_ss_n./(D0.^2);
    Cov(mask_final==0) = 0;
    figure('color', [1 1 1])
    subplot(3, 1, 1)
    imagesc(squeeze(Cov(:, :, slice)), [-sss sss]); axis image; axis off;
    colormap(mycmap)
    title('||Cov_ss||/D0^2')
    colorbar
    
    %% Powder average CTI fit
    %sel = [1 5 6 7 8];
    %gtab.bval1 = b1_all(sel);
    %gtab.bval2 = b2_all(sel);
    %gtab.cos2 = cos2ang(sel);
    
    gtab.bval1 = b1_all';
    gtab.bval2 = b2_all';
    gtab.cos2 = cos2ang';
    
    [DTOTAL_pa, KTOTAL_pa, KANISO_pa, KISO_pa, KINTRA_pa, S0_pa] = fun_CTIpa_rnh(Data_pa, gtab, ...
        mask_final, 0, 0, 0);
    
    [DTOTAL_mgc, KTOTAL_mgc, KANISO_mgc, KISO_mgc, S0_mgc] = fun_DIVIDE_rnh(Data_pa, gtab, ...
        mask_final, 0);
    
    mycmap = [0, 0, 0; parula];
    
    figure('color', [1 1 1])
    subplot(4, 1, 1)
    imagesc(squeeze(KTOTAL_pa(:, :, slice)), [0 1.5]); axis image; axis off;
    colormap(mycmap)
    colorbar
    
    subplot(4, 1, 2)
    imagesc(squeeze(KANISO_pa(:, :, slice)), [0 1.5]); axis image; axis off;
    colormap(mycmap)
    colorbar
    
    subplot(4, 1, 3)
    imagesc(squeeze(KISO_pa(:, :, slice)), [0 1.5]); axis image; axis off;
    colormap(mycmap)
    colorbar
    
    subplot(4, 1, 4)
    imagesc(squeeze(KINTRA_pa(:, :,slice)), [0 1.5]); axis image; axis off;
    colormap(mycmap)
    colorbar
    
    figure('color', [1 1 1])
    subplot(4, 1, 1)
    imagesc(squeeze(KANISO_mgc(:, :)+KISO_mgc(:, :)), [0 1.5]); axis image; axis off;
    colormap(mycmap)
    colorbar
    
    subplot(4, 1, 2)
    imagesc(squeeze(KANISO_mgc(:, :)), [0 1]); axis image; axis off;
    colormap(mycmap)
    colorbar
    
    subplot(4, 1, 3)
    imagesc(squeeze(KISO_mgc(:, :)), [0 1]); axis image; axis off;
    colormap(mycmap)
    colorbar
    
    
    if whichd == 1
        save_name = [dinfo(di).savename(1:end-7), 'cti', den, gst, cst];
    else
        save_name = [dinfo(di).savename(1:end-6), 'cti', den, gst, cst];
    end
    disp(save_name)
    save(save_name, ...
        'KISO', 'KANISO', 'KINTRA', 'KTOTAL', 'FA', 'MD', 'AD', 'RD', ...
        'MK', 'AK', 'RK', 'MKi', 'AKi', 'RKi', 'MKv', 'RKv', 'AKv',...
        'DTOTAL_pa', 'KTOTAL_pa', 'KANISO_pa', 'KISO_pa', 'KINTRA_pa', 'S0_pa',...
        'MKT', 'AKT', 'RKT', 'MKTv', 'AKTv', 'RKTv', 'MKTi', 'AKTi', 'RKTi', ...
        'DTOTAL_mgc', 'KTOTAL_mgc', 'KANISO_mgc', 'KISO_mgc', 'S0_mgc', ...
        'D0', 'D2', 'vW0', 'vW2', 'vW4', 'uW0', 'uW2', 'uW4', 'A0', 'A2',...
        'At', 'Cov_ss_n', 'Cov')
    
    
    close all
end
