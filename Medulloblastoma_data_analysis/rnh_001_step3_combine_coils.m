close all
clc

%addpath('..\DWIu_v1.11\DDE_toolbox')
%addpath('..\DWIu_v1.11\denoise')

ncoils = 4;

%% Section of code to edit
for whichd = 6:11 % Which data
    disp(whichd)
    
    %% Load file paths according to "which data" variable
    dinfo = fun_data_dirs(whichd);
    
    % check noise maps
    for di = 1:length(dinfo)
        for ci = 1:ncoils
            load([dinfo(di).savename, '_complex_den_', num2str(ci)], 'den_coil', 'gtab');
            
            if ci == 1
                SIZ = size(den_coil);
                data = abs(den_coil).^2;
            else
                data = data + abs(den_coil).^2;
            end
        end
        data = sqrt(data);
        save([dinfo(di).savename, '_complex_den'], 'data', 'gtab')
    end
end