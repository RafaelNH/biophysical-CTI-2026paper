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

denoised = true;
complex = true;
constraints = true;
Smooth_Gaussian = true;
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
if constraints
    cons = '_c';
else
    cons='';
end
if Smooth_Gaussian
    g = '_g';
else
g = '';
end

for whichd = 3:11
    dinfo = fun_data_dirs(whichd);
    if dinfo(1).folder>9
        load([dinfo(1).savename(1:end-6), 'MASKS.mat'], 'mask_final')
    else
        load([dinfo(1).savename(1:end-5), 'MASKS.mat'], 'mask_final')
    end
    
    for di = 1:length(dinfo)
        file_name = dinfo(di).savename;
        pt=[file_name, c, den, '_align_drift'];
        load(pt)
        
        if di == 1
            bvals1_all = gtab.bval1/1000;
            bvals2_all = gtab.bval2/1000;
            dir1_all = gtab.dir1;
            dir2_all = gtab.dir2;
            Data_all = data;
        else
            bvals1_new = gtab.bval1/1000;
            bvals2_new = gtab.bval2/1000;
            dir1_new = gtab.dir1;
            dir2_new = gtab.dir2;
            
            Data_all = cat(4, Data_all, data);
            bvals1_all = [bvals1_all, bvals1_new];
            bvals2_all = [bvals2_all, bvals2_new];
            dir1_all = [dir1_all, dir1_new];
            dir2_all = [dir2_all, dir2_new];
        end
    end
    
    
    dd = diag(dir1_all'*dir2_all);
    ddd = dd((bvals1_all + bvals2_all)> 0);
    
    figure, plot(bvals1_all), hold on, plot(bvals2_all, '--'), plot(dd - 2)
    
    if Smooth_Gaussian
        fwhm = 1.25;
        sd = fwhm/(sqrt(8*log(2))); % Convert FWHM to sd
        S = size(Data_all);
        Nvol = S(4);
        size_k = [7 7 7];
        for im = Nvol:-1:1
            for si = S(3):-1:1
                A = squeeze(Data_all(:,:,si,im));
                W = imgaussfilt(A, sd);
                Data_all(:,:,si,im) = W;
            end
        end
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
    %scal = (length(parula):-1:1)/length(parula);
    %mycmap = [scal', zeros(length(parula), 2); 0, 0, 0; parula];
    
    mycmap = turbo;
    for slice = [3, 5]
        %mycmap = [0, 0, 0; parula];
        %mycmap = parula;
        figure('color', [1 1 1])%, 'position', [-1400 -229 1000 840])
        subplot(1, 4, 1)
        imagesc(squeeze(KTOTAL(end:-1:1, :, slice)), [0 2]); axis image; axis off;
        colormap(mycmap)
        title('MK')
        colorbar
        
        subplot(1, 4, 2)
        imagesc(squeeze(KANISO(end:-1:1, :, slice)),  [0 1.5]); axis image; axis off;
        colormap(mycmap)
        title('Kaniso')
        colorbar
        
        subplot(1, 4, 3)
        imagesc(squeeze(KISO(end:-1:1, :, slice)),  [0 1.5]); axis image; axis off;
        colormap(mycmap)
        title('Kiso')
        colorbar
        
        subplot(1, 4, 4)
        imagesc(squeeze(KINTRA(end:-1:1, :, slice)), [0 1.5]); axis image; axis off;
        colormap(mycmap)
        title('Kintra')
        colorbar
    end
    
    [MD, AD, RD, FA, V1] = dti_metrics(pars0_all(:, :, :, 2:7), mask_final);
    
    MDs = MD(:, :, :);
    ADs = AD(:, :, :);
    RDs = RD(:, :, :);
    FAs = FA(:, :, :);
    for slice=3:7
        figure('color', [1 1 1], 'position', [-1400 -229 1000 840])
        subplot(4, 1, 1)
        imagesc(squeeze(MDs(end:-1:1, :, slice)), [0 1]); axis image; axis off;
        title('MD')
        colorbar
        
        subplot(4, 1, 2)
        imagesc(squeeze(ADs(end:-1:1 ,:, slice)), [0 1]); axis image; axis off;
        title('AD')
        colorbar
        
        subplot(4, 1, 3)
        imagesc(squeeze(RDs(end:-1:1, :, slice)), [0 1]); axis image; axis off;
        title('RD')
        colorbar
        
        subplot(4, 1, 4)
        imagesc(squeeze(FAs(end:-1:1, :, slice)), [0 1]); axis image; axis off;
        title('FA')
        colorbar
        
        colormap(gray)
    end
    
    
    
    [DT, KT, KTi, KTv] = fun_resolve_kurtosis_tensors(pars0_all);
    
    [MKT, AKT, RKT] = dki_tensor_metrics(DT, KT, mask_final);
    [MKTv, AKTv, RKTv] = dki_tensor_metrics(DT, KTv, mask_final);
    [MKTi, AKTi, RKTi] = dki_tensor_metrics(DT, KTi, mask_final);
    
    
    sss = 1.5;
    
    for slice=[3,5]
        figure('color', [1 1 1])%, 'position', [-1400 -229 1000 840])
        subplot(3, 3, 3)
        imagesc(squeeze(RKT(end:-1:1, :, slice)), [0 sss]); axis image; axis off;
        title('RKT')
        
        subplot(3, 3, 2)
        imagesc(squeeze(AKT(end:-1:1, :, slice)), [0 sss]); axis image; axis off;
        title('AKT')
        
        subplot(3, 3, 1)
        imagesc(squeeze(MKT(end:-1:1, :, slice)), [0 sss]); axis image; axis off;
        title('MKT')
        
        
        subplot(3, 3, 6)
        imagesc(squeeze(RKTv(end:-1:1, :,  slice)), [0 sss]); axis image; axis off;
        title('RKTv')
        
        subplot(3, 3, 5)
        imagesc(squeeze(AKTv(end:-1:1, :, slice)), [0 sss]); axis image; axis off;
        title('AKTv')
        
        subplot(3, 3, 4)
        imagesc(squeeze(MKTv(end:-1:1, :, slice)), [0 sss]); axis image; axis off;
        title('MKTv')
        
        
        subplot(3, 3, 9)
        imagesc(squeeze(RKTi(end:-1:1, :, slice)), [0 sss]); axis image; axis off;
        title('RKTi')
        
        subplot(3, 3, 8)
        imagesc(squeeze(AKTi(end:-1:1, :, slice)), [0 sss]); axis image; axis off;
        title('AKTi')
        
        subplot(3, 3, 7)
        imagesc(squeeze(MKTi(end:-1:1, :, slice)), [0 sss]); axis image; axis off;
        title('MKTi')
        colormap(mycmap)
    end
    
    %% RICE
    %[D2mat, S0mat, S2mat, S4mat, A0mat, A2mat, SSV, COVi] = fun_resolve_invariants(pars0_all, mask_final);
    [D0, D2, vW0, vW2, vW4, uW0, uW2, uW4, A0, A2, At, Cov_ss_n] = fun_resolve_cti_invariants(pars0_all, mask_final);
    
    for slice=6:7
        Cov = Cov_ss_n;
        Cov(mask_final==0) = 0;
        figure('color', [1 1 1])
        subplot(1, 1, 1)
        imagesc(squeeze(sqrt(Cov(end:-1:1, :, slice))), [0 0.6]); axis image; axis off;
        colormap(mycmap)
        title('||Cov_ss||')
        colorbar
    end
    
    file_name = dinfo(1).savename;
    if dinfo(1).folder>9
        save([dinfo(1).savename(1:end-6), 'cti_', c, den, 'align_drift', g, cons], ...
            'KISO', 'KANISO', 'KINTRA', 'KTOTAL', 'FA', 'MD', 'AD', 'RD', 'S0',... ...
            'MKT', 'AKT', 'RKT', 'MKTv', 'AKTv', 'RKTv', 'MKTi', 'AKTi', 'RKTi', ...
            'D0', 'D2', 'vW0', 'vW2', 'vW4', 'uW0', 'uW2', 'uW4', 'A0', 'A2',...
            'At', 'Cov_ss_n')
    else
        save([dinfo(1).savename(1:end-5), 'cti_', c, den, 'align_drift', g, cons], ...
            'KISO', 'KANISO', 'KINTRA', 'KTOTAL', 'FA', 'MD', 'AD', 'RD', 'S0',... ...
            'MKT', 'AKT', 'RKT', 'MKTv', 'AKTv', 'RKTv', 'MKTi', 'AKTi', 'RKTi', ...
            'D0', 'D2', 'vW0', 'vW2', 'vW4', 'uW0', 'uW2', 'uW4', 'A0', 'A2',...
            'At', 'Cov_ss_n')
    end
    
end
