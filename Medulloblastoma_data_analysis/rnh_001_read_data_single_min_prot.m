%% Prepare CTI data acquired with a minimal protocol using a single
%% experiment (note this assumed hard coded order of acquisitions).

% Notes: I am excluding the first b0

clear all
close all
clc

addpath('..\DWIu_v1.11\DDE_toolbox')
addpath(genpath('..\DWIu_v1.11\pvtools\'));
addpath(genpath('..\DWIu_v1.11\masktools\'));

%% Section of code to edit
whichd = 11; % Which data
dinfo = fun_data_dirs(whichd);

read_complex = true;
if read_complex
    complex = '_complex';
else
    complex = '';
end

for di = 1:length(dinfo)
    
    %% Read Bruker files using pvtools
    
    file_dir = [dinfo(di).name, filesep, 'pdata'];
    
    if read_complex
        directories = dir(file_dir);
        proc_num = str2num(directories(end).name);
    else
        proc_num = 1;
    end
    
    source = strcat(file_dir, filesep, num2str(proc_num), filesep);
    
    visu_params = readBrukerParamFile(strcat(source, 'visu_pars'));
    
    acq_params = readBrukerParamFile(strcat(dinfo(di).name, filesep,'acqp'));
    
    method_params = readBrukerParamFile(strcat(dinfo(di).name, filesep, 'method'));
    
    disp(method_params.RNH_Angle)
    disp(method_params.PVM_NAverages)
    % load data depending on the number of tmix acquired
    
    [data_all, Visu] = readBruker2dseq(strcat(source,'2dseq'), visu_params);
    
    permuteOrder=1:ndims(data_all);
    permuteOrder([1 2])=permuteOrder([2 1]);
    data_all = flip(permute(squeeze(data_all),permuteOrder),2);  %% sometimes needed
    if read_complex
        ncoils = size(data_all, length(size(data_all)));
        nvol_all = size(data_all, length(size(data_all))-1);
    else
        nvol_all = size(data_all, length(size(data_all)));
    end
    
    %% Read acquisition parameters
    bval_a = method_params.AI_DWTotalBvaluea;
    bval_b = method_params.AI_DWTotalBvalueb;
    bval_c = method_params.AI_DWTotalBvaluec;
    bval_d = method_params.AI_DWTotalBvalued;
    nb0 = method_params.AI_NB0;
    rang = method_params.RNH_Angle;
    
    % Read gradient time information
    sdelta = method_params.AI_DWGradDur1*0.001;
    bdelta = method_params.AI_DWGradSep1;
    tmix = method_params.AI_DWGradMix;
    nvol = nvol_all/length(tmix);
    
    % I am excluding the first nb0
    
    % Read gradient directions
    dir2 = method_params.AI_DWGradDir2(2:nvol, :);
    dir1 = method_params.AI_DWGradDir1(2:nvol, :);
    
    if read_complex
        data_all = data_all(:, :, :, 2:nvol, :);
    else
        data_all = data_all(:, :, :, 2:nvol);
    end
    g2a = method_params.AI_DWGradStr2a;
    g2b = method_params.AI_DWGradStr2b;
    g2c = method_params.AI_DWGradStr2c;
    g2d = method_params.AI_DWGradStr2d;
    g1a = method_params.AI_DWGradStr1a;
    g1b = method_params.AI_DWGradStr1b;
    g1c = method_params.AI_DWGradStr1c;
    g1d = method_params.AI_DWGradStr1d;
    
    % convert gradient intensities to bvalues
    const_a = bval_a / (g1a * g1a + g2a * g2a);
    const_b = bval_b / (g1b * g1b + g2b * g2b);
    const_c = bval_c / (g1c * g1c + g2c * g2c);
    const_d = bval_d / (g1d * g1d + g2d * g2d);
    
    b1_a = const_a * g1a * g1a;
    b2_a = const_a * g2a * g2a;
    
    b1_b = const_b * g1b * g1b;
    b2_b = const_b * g2b * g2b;
    
    b1_c = const_c * g1c * g1c;
    b2_c = const_c * g2c * g2c;
    
    b1_d = const_d * g1d * g1d;
    b2_d = const_d * g2d * g2d;
    
    % Plot b-values to inspect acquired parameters
    disp('Acquired b-values (b1, b2, bt, bt_check):')
    disp([b1_a, b2_a, b1_a + b2_a, bval_a;
        b1_b, b2_b, b1_b + b2_b, bval_b;
        b1_c, b2_c, b1_c + b2_c, bval_c;
        b1_d, b2_d, b1_d + b2_d, bval_d]);
    
    % calculate summary parameters
    b1_all = [b1_a, b1_b, b1_c, b1_d];
    b2_all = [b2_a, b2_b, b2_c, b2_d];
    
    % Reconstruct gtab
    gtab.bval1 = method_params.RNH_bvals1(2:nvol);
    gtab.bval2 = method_params.RNH_bvals2(2:nvol);
    gtab.bval = gtab.bval1 + gtab.bval2;
    gtab.dir1 = dir1';  % all directions 1
    gtab.dir2 = dir2';  % all directions 2
    gtab.rang = rang;
    gtab.sdelta = sdelta;
    gtab.bdelta = bdelta;
    gtab.tmix = tmix;
    gtab.uc2a = [1, 1, 1, cos(rang/180*pi).^2, 1];
    gtab.ub1 = [b2_a, b1_a, b1_b, b1_c, b1_d];
    gtab.ub2 = [b2_a, b2_a, b2_b, b2_c, b2_d];
    gtab.ub = gtab.ub1 + gtab.ub2;
    
    % shells information used for data correction
    gtab.shells = gtab.bval*0;
    gtab.shells(round(gtab.bval1)> round(gtab.bval2)) = 2;
    gtab.shells(round(gtab.bval1) == round(gtab.bval2)) = 3;
    gtab.shells(round(diag(dir1*dir2'))==0)=4;
    gtab.shells(round(gtab.bval) == round(b1_d + b2_d)) = 5;
    gtab.shells(round(gtab.bval)==round(b2_a*2)) = 1;
    
    
    figure
    plot(gtab.shells)
    
    
    if read_complex
        for ci = 1:ncoils
            data = squeeze(data_all(:, :, :, :, ci));
            save([dinfo(di).savename, complex, '_', num2str(ci)] , 'data', 'gtab');
        end
    else
        data = data_all;
        save(dinfo(di).savename , 'data', 'gtab');
    end
end
