%% Prepare CTI data acquired with a minimal protocol using a single
%% experiment (note this assumed hard coded order of acquisitions).

close all
clc

addpath('..\DWIu_toolbox_v1.11\DDE_toolbox')
addpath(genpath('..\DWIu_toolbox_v1.11\pvtools\'));

%% Section of code to edit
for whichd = 10; % Which data
proc_num = 1; % specify which processing folder you are interested in

%% Load file paths according to "which data" variable
dinfo = fun_data_dirs(whichd);

for di = 1:length(dinfo)

%% Read Bruker files using pvtools

source = strcat(dinfo(di).name, filesep, ...
    'pdata', filesep, num2str(proc_num), filesep);

visu_params = readBrukerParamFile(strcat(source, 'visu_pars'));

acq_params = readBrukerParamFile(strcat(dinfo(di).name, filesep,'acqp'));

method_params = readBrukerParamFile(strcat(dinfo(di).name, filesep, 'method'));
method_params.Method
[imagee, Visu] = readBruker2dseq(strcat(source,'2dseq'), visu_params);

% Orientation may need
permuteOrder=1:ndims(imagee);
permuteOrder([1 2])=permuteOrder([2 1]);
data = flip(permute(squeeze(imagee),permuteOrder),2);  %% sometimes needed
nvol = size(data, 4);

figure('color', [1 1 1], 'position', [89, 79, 1100, 620])
if di == 1
    subplot(1, 3, 1)
    S0 = squeeze(data(:, :, 10, 1));
    imagesc(S0);
    S0ref = S0;
    axis image
    axis off
else
    subplot(1, 3, 1)
    S0 = squeeze(data(:, :, 10, 1));
    imagesc(S0);
    axis image
    axis off
    subplot(1, 3, 2)
    imagesc(S0ref);
    axis image
    axis off
    subplot(1, 3, 3)
    imagesc(S0 - S0ref, [-5 5]);
    axis image
    axis off
end


%% Read acquisition parameters
bval = method_params.AI_DWTotalBvalue;
nb0 = method_params.AI_NB0;
rang = method_params.RNH_Angle;

% Read gradient directions
dir2 = method_params.AI_DWGradDir2(1:nvol, :);
dir1 = method_params.AI_DWGradDir1(1:nvol, :);

% Read gradient information
sdelta = method_params.AI_DWGradDur1*0.001;
bdelta = method_params.AI_DWGradSep1;
tmix = method_params.AI_DWGradMix;

g1 = method_params.AI_DWGradStr1;
g2 = method_params.AI_DWGradStr2;


% convert gradient intensities to bvalues
const = bval / (g1 * g1 + g2 * g2);

b1 = const * g1 * g1;
b2 = const * g2 * g2;

% Plot b-values to inspect acquired parameters
disp('Acquired b-values (b1, b2, bt, bt_check):')
disp([b1, b2, b1 + b2, bval]);

% calculate summary parameters
b1_all = b1;
b2_all = b2;
gtab = fun_create_DDE_gtab_folder(b1_all, b2_all, dir1, dir2, nb0, rang, ...
    sdelta, bdelta, tmix);

disp('tmix:')
disp(tmix)
disp('rang:')
disp(rang)

drawnow
pause(0.05)
%% save data
save(dinfo(di).savename, 'data', 'gtab');
end
disp(di)


% some code for data inspection
% for boi = 1:20
%     figure
%     subplot(1, 3, 1)
%     S0 = squeeze(data(:, :, 10, boi));
%     S0_ref = squeeze(data_ref(:, :, 10, boi));
%     imagesc(S0);
%     axis image
%     axis off
%     subplot(1, 3, 2)
%     imagesc(S0_ref);
%     axis image
%     axis off
%     subplot(1, 3, 3)
%     imagesc(S0 - S0_ref, [-5 5]);
%     axis image
%     axis off
% end
% % 
% % 
% for si = 1:20
%     figure
%     subplot(1, 3, 1)
%     S0 = squeeze(data(:, :, si, 1));
%     S0_ref = squeeze(data_ref(:, :, si, 1));
%     imagesc(S0);
%     axis image
%     axis off
%     subplot(1, 3, 2)
%     imagesc(S0_ref);
%     axis image
%     axis off
%     subplot(1, 3, 3)
%     imagesc(S0 - S0_ref, [-5 5]);
%     axis image
%     axis off
% end
end