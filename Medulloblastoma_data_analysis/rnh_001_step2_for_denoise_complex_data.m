close all
clc

addpath('..\DWIu_v1.11\DDE_toolbox')
addpath('..\DWIu_v1.11\denoise')

ncoils = 4;

%% Section of code to edit
for whichd = 8:11 % Which data
    
    %% Load file paths according to "which data" variable
    dinfo = fun_data_dirs(whichd);
    
    % check noise maps
    for di = 1:length(dinfo)
        for ci = 1:ncoils
            disp([dinfo(di).savename,'_complex_', num2str(ci)])
            load([dinfo(di).savename,'_complex_', num2str(ci)], 'data', 'gtab');
            
            S0_coil = data(:, :, :, gtab.shells==1);
            S2_coil = var(S0_coil(:, :, :, end-10:end), 0, 4);
            mS0_coil = mean(abs(S0_coil), 4);
            
            figure('color', [1 1 1])
            subplot(4, 1, 1)
            imagesc(mS0_coil(:, :), [0, 4]);
            subplot(4, 1, 2)
            imagesc(S2_coil(:, :), [0 1]);

            [den_coil, S2new, P] = denoise_tpca(data, [11 11 3], S2_coil);
            
            %den_complex
            figure('color', [1 1 1])
            
            subplot(3, 2, 1)
            data_sel = abs(data(end:-1:1, :, 6:9, 21));
            imagesc(data_sel(:, :), [0 20]);
            
            subplot(3, 2, 3)
            den_sel = abs(den_coil(end:-1:1, :, 6:9, 21));
            imagesc(den_sel(:, :), [0 20]);
            
            subplot(3, 2, 4)
            imagesc(data_sel(:, :)-den_sel(:,:), [-2 2]);
            
            subplot(3, 2, 5)
            P_sel = P(:, :, 6:9);
            imagesc(P_sel(:, :), [0 30]);
            colorbar
            
            subplot(3, 2, 6)
            S2n_sel = S2new(:, :, 6:9);
            imagesc(S2n_sel(:, :), [0 5]);
            drawnow
            
            save([dinfo(di).savename, '_complex_den_', num2str(ci)], 'den_coil', 'gtab', 'P', 'S2new');
        end
    end
end