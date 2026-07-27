clear all
close all
clc

fs = filesep;
addpath(['..', fs, 'DWIu_toolbox_v1.11', fs, 'Spherical_functions_and_directions/'])
addpath(['..',fs,'NIFTI_toolbox'])
addpath(['..',fs,'DWIu_toolbox_v1.11',fs,'DDE_toolbox'])
addpath(['..',fs,'DWIu_toolbox_v1.11',fs,'Spherical_functions_and_directions'])
load('dirs45.mat')

Smooth_Gaussian = true;
const = true;

for whichd = 10
    [path, dwi_name, bval_name, bvec_name, mask_name, datatype] = fun_data_dir(whichd);
    
    % Load data
    v = load_untouch_nii([path, dwi_name]);
    data = double(v.img);
    
    v = load_untouch_nii([path, mask_name]);
    mask_final = double(v.img);
    %mask_final(:, :, 1:33) = 0;
    %mask_final(:, :, 35:end) = 0;
    
    
    if Smooth_Gaussian
        fwhm = 1.25;
        sd = fwhm/(sqrt(8*log(2))); % Convert FWHM to sd
        S = size(data);
        Nvol = S(4);
        size_k = [7 7 7];
        for im = Nvol:-1:1
            A = squeeze(data(:,:,:,im));
            W = smooth3(A,'gaussian',size_k,sd);
            data(:,:,:,im) = W;
        end
    end
    
    data(data(:)<0) = 0;
    
    bvals = load([path, bval_name]);
    bvecs = load([path, bvec_name]);
    
    
    gtab = fun_reconst_gtab_for_human(bvals/1000, bvecs, datatype);
    
    if whichd == 10
    save('gtab', 'gtab')
    end
    
    S0s = data(:, :, :, gtab.bval==0);
    S0 = mean(S0s, 4);
    
    figure('color', [1 1 1])
    % for nv = 1:45
    %     for s = 1:60
    %     imagesc(squeeze(data(:, :, s, nv)), [-10 0])
    %     title(nv)
    %     drawnow
    %     end
    % end
    % figure('color', [1 1 1])
    % for s = 1:60
    %     imagesc(squeeze(S0(:, :, s)))
    %     title(nv)
    %     drawnow
    % end
    %subplot(2, 1, 1)
    %imagesc(squeeze(S0s(:, :, 30, 1)))
    %subplot(2, 1, 2)
    %imagesc(mask_final(:, :))
    
    [Data_pa] = fun_compute_powderaverage(data, gtab);
    
    bvals1_all = gtab.bval1;
    bvals2_all = gtab.bval2;
    dir1_all = gtab.dir1;
    dir2_all = gtab.dir2;
    b1_all = [0, gtab.ub1];
    b2_all = [0, gtab.ub2];
    cos2ang = [1, gtab.uc2a];
    
    if const
        %         [D11,D22,D33,D12,D13,D23,...
        %             W1111, W2222, W3333, W1112, W1113,...
        %             W1222, W2223, W1333, W2333, W1122,...
        %             W1133, W2233, W1123, W1223, W1233,...
        %             K1111, K2222, K3333, K1112, K1113, ...
        %             K1222, K2223, K1333, K2333, K1122, ...
        %             K1133, K2233, K1123, K1322, K1233, ...
        %             K1212, K1313, K2323, K1223, K1323, K1213, S0, SSE, Ndiff, Nkurt, Ncov] = ...
        %             fun_DDE_CLLS_comp_b2_rh(data, mask_final, ...
        %             bvals1_all', bvals2_all', dir1_all, dir2_all);
        [D11,D22,D33,D12,D13,D23,...
            W1111, W2222, W3333, W1112, W1113,...
            W1222, W2223, W1333, W2333, W1122,...
            W1133, W2233, W1123, W1223, W1233,...
            K1111, K2222, K3333, K1112, K1113, ...
            K1222, K2223, K1333, K2333, K1122, ...
            K1133, K2233, K1123, K1322, K1233, ...
            K1212, K1313, K2323, K1223, K1323, K1213, S0, SSE, Ndiff, Nkurt, Ncov] = ...
            fun_DDE_CLLS_comp_b2_improved_rh(data, mask_final, ...
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
            fun_DDE_ULLS_comp_b2_rh(data, mask_final, ...
            bvals1_all', bvals2_all', dir1_all, dir2_all);
    end
    
    N = size(data);
    
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
    slice = 34;
    mycmap = [0, 0, 0; parula];
    
    figure('color', [1 1 1])
    subplot(1, 4, 1)
    imagesc(squeeze(KTOTAL(:, end:-1:1, slice))', [0 1.5]); axis image; axis off;
    colormap(mycmap)
    colorbar
    
    subplot(1, 4, 2)
    imagesc(squeeze(KANISO(:, end:-1:1, slice))', [0 1]); axis image; axis off;
    colormap(mycmap)
    colorbar
    
    subplot(1, 4, 3)
    imagesc(squeeze(KISO(:, end:-1:1, slice))', [0 1]); axis image; axis off;
    colormap(mycmap)
    colorbar
    
    subplot(1, 4, 4)
    imagesc(squeeze(KINTRA(:, end:-1:1, slice))', [0 1]); axis image; axis off;
    colormap(mycmap)
    colorbar
    
    %% Diffusion Tensor Metrics
    [MD, AD, RD, FA, V1] = dti_metrics(pars0_all(:, :, :, 2:7), mask_final);
    
    MDs = MD(:, :, :);
    ADs = AD(:, :, :);
    RDs = RD(:, :, :);
    FAs = FA(:, :, :);
    
    figure('color', [1 1 1])
    subplot(1, 4, 1)
    imagesc(squeeze(MDs(:, end:-1:1, slice))', [0 2]); axis image; axis off;
    title('MD')
    colorbar
    
    subplot(1, 4, 2)
    imagesc(squeeze(ADs(:, end:-1:1, slice))', [0 2]); axis image; axis off;
    title('AD')
    colorbar
    
    subplot(1, 4, 3)
    imagesc(squeeze(RDs(:, end:-1:1, slice))', [0 2]); axis image; axis off;
    title('RD')
    colorbar
    
    subplot(1, 4, 4)
    imagesc(squeeze(FAs(:, end:-1:1, slice))', [0 0.7]); axis image; axis off;
    title('FA')
    colorbar
    
    colormap(gray)
    
    %% Resolve micro and variance kurtosis tensors
    [DT, KT, KTi, KTv] = fun_resolve_kurtosis_tensors(pars0_all);
    
    
    % Kurtosis tensor metrics
    [MKT, AKT, RKT] = dki_tensor_metrics(DT, KT, mask_final);
    % Variance kurtosis tensor metrics
    [MKTv, AKTv, RKTv] = dki_tensor_metrics(DT, KTv, mask_final);
    % Micro kurtosis tensor metrics
    [MKTi, AKTi, RKTi] = dki_tensor_metrics(DT, KTi, mask_final);
    
    slide = 34;
    
    %% Powder average CTI fit
    
    gtab.bval1 = b1_all';
    gtab.bval2 = b2_all';
    gtab.cos2 = cos2ang';
    
    [DTOTAL_pa, KTOTAL_pa, KANISO_pa, KISO_pa, KINTRA_pa, S0_pa] = fun_CTIpa_rnh(Data_pa, gtab, ...
        mask_final, 0, 0, 0);
    
    mycmap = [0, 0, 0; parula];
    
    figure('color', [1 1 1])
    subplot(1, 4, 1)
    imagesc(squeeze(KTOTAL_pa(:, end:-1:1, slice))', [0 1.5]); axis image; axis off;
    colormap(mycmap)
    colorbar
    
    subplot(1, 4, 2)
    imagesc(squeeze(KANISO_pa(:, end:-1:1, slice))', [0 1]); axis image; axis off;
    colormap(mycmap)
    colorbar
    
    subplot(1, 4, 3)
    imagesc(squeeze(KISO_pa(:, end:-1:1, slice))', [0 1]); axis image; axis off;
    colormap(mycmap)
    colorbar
    
    subplot(1, 4, 4)
    imagesc(squeeze(KINTRA_pa(:, end:-1:1, slice))', [0 1]); axis image; axis off;
    colormap(mycmap)
    colorbar
    
    if Smooth_Gaussian
        if const
            save([path, fs, 'cti_gs_c2'], ...
                'KISO', 'KANISO', 'KINTRA', 'KTOTAL', 'FA', 'MD', 'AD', 'RD', 'V1', ...
                'DTOTAL_pa', 'KTOTAL_pa', 'KANISO_pa', 'KISO_pa', 'KINTRA_pa', 'S0_pa',...
                'MKT', 'AKT', 'RKT', 'MKTv', 'AKTv', 'RKTv', 'MKTi', 'AKTi', 'RKTi', 'S0', ...
                'Ndiff', 'Nkurt', 'Ncov');
            
            path_metrics = [path, fs, 'CTImetrics_gs_c2'];
            
        else
            save([path, fs, 'cti_gs'], ...
                'KISO', 'KANISO', 'KINTRA', 'KTOTAL', 'FA', 'MD', 'AD', 'RD', 'V1', ...
                'DTOTAL_pa', 'KTOTAL_pa', 'KANISO_pa', 'KISO_pa', 'KINTRA_pa', 'S0_pa',...
                'MKT', 'AKT', 'RKT', 'MKTv', 'AKTv', 'RKTv', 'MKTi', 'AKTi', 'RKTi', 'S0');
            
            path_metrics = [path, fs, 'CTImetrics_gs'];
        end
    else
        if const
            save([path, fs, 'cti_c'], ...
                'KISO', 'KANISO', 'KINTRA', 'KTOTAL', 'FA', 'MD', 'AD', 'RD', 'V1', ...
                'DTOTAL_pa', 'KTOTAL_pa', 'KANISO_pa', 'KISO_pa', 'KINTRA_pa', 'S0_pa',...
                'MKT', 'AKT', 'RKT', 'MKTv', 'AKTv', 'RKTv', 'MKTi', 'AKTi', 'RKTi', 'S0', 'Ndiff', 'Nkurt', 'Ncov');
            
            path_metrics = [path, fs, 'CTImetrics_c'];
            
        else
            save([path, fs, 'cti'], ...
                'KISO', 'KANISO', 'KINTRA', 'KTOTAL', 'FA', 'MD', 'AD', 'RD', 'V1', ...
                'DTOTAL_pa', 'KTOTAL_pa', 'KANISO_pa', 'KISO_pa', 'KINTRA_pa', 'S0_pa',...
                'MKT', 'AKT', 'RKT', 'MKTv', 'AKTv', 'RKTv', 'MKTi', 'AKTi', 'RKTi', 'S0');
            
            path_metrics = [path, fs, 'CTImetrics'];
        end
    end
    
    mkdir(path_metrics)
    vnii = v;
    fun_save_nii(S0, vnii, [path_metrics, fs, 'S0.nii'])
    fun_save_nii(KISO, vnii, [path_metrics, fs, 'KISO.nii'])
    fun_save_nii(KANISO, vnii, [path_metrics, fs, 'KANISO.nii'])
    fun_save_nii(KINTRA, vnii, [path_metrics, fs, 'KINTRA.nii'])
    fun_save_nii(KTOTAL, vnii, [path_metrics, fs, 'KTOTAL.nii'])
    fun_save_nii(FA, vnii, [path_metrics, fs, 'FA.nii'])
    fun_save_nii(MD, vnii, [path_metrics, fs, 'MD.nii'])
    fun_save_nii(AD, vnii, [path_metrics, fs, 'AD.nii'])
    fun_save_nii(RD, vnii, [path_metrics, fs, 'RD.nii'])
    fun_save_nii(MKT, vnii, [path_metrics, fs, 'MKT.nii'])
    fun_save_nii(AKT, vnii, [path_metrics, fs, 'AKT.nii'])
    fun_save_nii(RKT, vnii, [path_metrics, fs, 'RKT.nii'])
    fun_save_nii(MKTv, vnii, [path_metrics, fs, 'MKTv.nii'])
    fun_save_nii(AKTv, vnii, [path_metrics, fs, 'AKTv.nii'])
    fun_save_nii(RKTv, vnii, [path_metrics, fs, 'RKTv.nii'])
    fun_save_nii(MKTi, vnii, [path_metrics, fs, 'MKTi.nii'])
    fun_save_nii(AKTi, vnii, [path_metrics, fs, 'AKTi.nii'])
    fun_save_nii(RKTi, vnii, [path_metrics, fs, 'RKTi.nii'])
    %fun_save_nii(Cov_ss_n, vnii, [path_metrics, fs, 'Cov_ss_n.nii'])
    close all
end