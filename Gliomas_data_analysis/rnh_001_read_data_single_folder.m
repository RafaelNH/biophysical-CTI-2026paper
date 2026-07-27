%% Prepare CTI data acquired with a minimal protocol using a single
%% experiment (note this assumed hard coded order of acquisitions).

clear all
close all
clc

addpath('..\DWIu_v1.11\DDE_toolbox')
addpath(genpath('..\DWIu_v1.11\pvtools\'));
addpath(genpath('..\DWIu_v1.11\masktools\'));

%% Section of code to edit
whichd = 7; % Which data
dinfo = fun_data_dirs(whichd);
proc_num = 1; % specify which processing folder you are interested in

for di = 1:length(dinfo)
    
    %% Read Bruker files using pvtools
    
    source = strcat(dinfo(di).name, filesep, ...
        'pdata', filesep, num2str(proc_num), filesep);
    
    visu_params = readBrukerParamFile(strcat(source, 'visu_pars'));
    
    acq_params = readBrukerParamFile(strcat(dinfo(di).name, filesep,'acqp'));
    
    method_params = readBrukerParamFile(strcat(dinfo(di).name, filesep, 'method'));
    
    % load data depending on the number of tmix acquired
    
    [data_all, Visu] = readBruker2dseq(strcat(source,'2dseq'), visu_params);
    
    permuteOrder=1:ndims(data_all);
    permuteOrder([1 2])=permuteOrder([2 1]);
    data_all = flip(permute(squeeze(data_all),permuteOrder),2);  %% sometimes needed
    nvol_all = size(data_all, 4);
    
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
    
    % Read gradient directions
    dir2 = method_params.AI_DWGradDir2(1:nvol, :);
    dir1 = method_params.AI_DWGradDir1(1:nvol, :);
    
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
    
    if method_params.Method == "<User:ra_DDE_cti_all>"
        gtab = fun_create_DDE_gtab_expanded(b1_all, b2_all, dir1, dir2,...
            0, 90, sdelta, bdelta, tmix);
        
        data = data_all;
        save(dinfo(di).savename, 'data', 'gtab');
    elseif method_params.Method == "<User:rnh_DDE_cti_all>"
        gtab = fun_create_DDE_gtab_expanded(b1_all, b2_all, dir1, dir2, nb0, rang, ...
            sdelta, bdelta, tmix);
        
        data = data_all;
        save(dinfo(di).savename, 'data', 'gtab');
    end
end
