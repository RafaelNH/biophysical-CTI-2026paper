close all
clc

addpath('..\DWIu_v1.11\DDE_toolbox')
addpath('..\DWIu_v1.11\denoise')

%% Section of code to edit
for whichd = 2:8 % Which data
    proc_num = 1; % specify which processing folder you are interested in
    
    %% Load file paths according to "which data" variable
    dinfo = fun_data_dirs(whichd);
    
    % check noise maps
    for di = 1:length(dinfo)
        load(dinfo(di).savename, 'data', 'gtab');
        bvals = gtab.bval1 + gtab.bval2;
        S0 = data(:, :, :, bvals==0);
        S0 = S0(:, :, :, 1:10);
        S2 = var(S0, 0, 4);
    
        figure('color', [1 1 1])
        S2_sel = S2(:, :, 1:8);
        subplot(3, 1, 1)
        imagesc(S2_sel(:, :), [0 50]);
        S2_sel = S2(:, :, 9:16);
        subplot(3, 1, 2)
        imagesc(S2_sel(:, :), [0 50]);

    end
    
   close all
    
    % denoise
    for di = 1:length(dinfo)
        load(dinfo(di).savename, 'data', 'gtab');
        bvals = gtab.bval1 + gtab.bval2;
        S0 = data(:, :, :, bvals==0);
        S0 = S0(:, :, :, 1:10);
        S2 = var(S0, 0, 4);
        
        [den, S2new, P] = denoise_tpca(data, [11 11 6], S2);
        
%         figure('color', [1 1 1])
%         
%         subplot(3, 2, 1)
%         data_sel = data(:, :, 6:10, 21);
%         imagesc(data_sel(:, :), [0 20]);
%         
%         subplot(3, 2, 3)
%         den_sel = den(:, :, 6:10, 21);
%         imagesc(den_sel(:, :), [0 20]);
%         
%         subplot(3, 2, 4)
%         imagesc(data_sel(:, :)-den_sel(:,:), [-2 2]);
%         
%         subplot(3, 2, 5)
%         P_sel = P(:, :, 6:10);
%         imagesc(P_sel(:, :), [0 30]);
%         colorbar
%         
%         subplot(3, 2, 6)
%         S2n_sel = S2new(:, :, 6:10);
%         imagesc(S2n_sel(:, :), [0 5]);
%         drawnow
%         pause(0.01)
        
        data = den;
        save([dinfo(di).savename, '_den'], 'data', 'gtab', 'P', 'S2new');
    end
    disp(whichd)
end