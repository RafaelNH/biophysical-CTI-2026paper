close all
clear all
clc

fs = filesep;
addpath(['..', fs, 'DWIu_toolbox_v1.11', fs, 'DKI_toolbox/'])
addpath(['..', fs, 'DWIu_toolbox_v1.11', fs, 'DDE_toolbox/'])
addpath(['..', fs, 'DWIu_toolbox_v1.11', fs, 'Spherical_functions_and_directions/'])

As = 0;
as = As*100;
seed = 123;
intraD = 250;
extraD = 150;

% first aspect - compare diffusion between plates and diffusion inside cylinders
times = 0.1:0.1:10;
for ti = length(times):-1:1
    [dc(ti), kc(ti), wc(ti)] = diffusion_in_cylinder(times(ti), 9, 20, intraD/100);
    [dp(ti), kp(ti), wp(ti)] = diffusion_in_plane(times(ti), 18, 20, intraD/100);
end

figure
subplot(1, 2, 1)
plot(times, dc)
hold on
plot(times, dp)

subplot(1, 2, 2)
plot(times, kc)
hold on
plot(times, kp)


pth = ['MC_results_from_step001', fs];
load([pth,'tensors_as',num2str(As),'_intraD',num2str(intraD),'_seed',(num2str(seed))]);

% second aspect - check if we can reproduce the the results of the first
% simulation scenario using the analytical solution R=1

vf = 0.7;
load([pth,'tensors_as',num2str(as),'_f',num2str(vf*100),'_extraD',num2str(extraD),'_seed',(num2str(seed))]);
Evector = [0, 0, 1;...
    0, 1, 0;...
    1, 0, 0];
[dt_ex_unit, wt_ex_unit] = fun_rotate_tensors_DKI(dt_ex, wt_ex, Evector);
mde = mean(dt_ex_unit(1:3));
[dt_in_unit, wt_in_unit, mdi] = compute_cylinder_unit_tensors(10, 1, intraD/100, 20);
md = mdi * vf + mde*(1-vf);

load('dirs1024.mat')
disp('sampling bingham...')
ODI = 0.12;
k = 1/tan(pi/2*ODI);
angle = asin(sqrt(1/(2*k)))/pi*180;

V2 = V;
V3 = V;

fp = bingham(V, 1, k, k, [0, 1, 0]', [0, 0, 1]');

for vi = 1:size(V, 1)
    Xa = V(vi, :);
    tang = 2*pi*rand();
    pang = acos(1 - 2*rand());
    Xri = sin(pang) * cos(tang);
    Xrj = sin(pang) * sin(tang);
    Xrk = cos(pang);
    
    while abs(Xa(1) - Xri) < 0.05
        tang = 2*pi*rand();
        pang = acos(1 - 2*rand());
        Xri = sin(pang) * cos(tang);
        Xrj = sin(pang) * sin(tang);
        Xrk = cos(pang);
    end
    Xr = [Xri, Xrj, Xrk];
    
    Xb=cross(Xa,Xr);
    Xb=Xb/norm(Xb);
    Xc=cross(Xa,Xb);
    Xc=Xc/norm(Xc);
    
    V2(vi, :) = Xb;
    V3(vi, :) = Xc;
    
end

fp = fp / sum(fp);
fp_all = fp; % Since nbins = 1, fp_all is exactly equal to fp

% Process Intracellular Compartment
disp('Processing intracellular compartment...')
[uwt_in, dt_in_total, dt_in_all] = compute_compartment_tensors(V, V2, V3, fp, dt_in_unit, wt_in_unit);

% Process Extracellular Compartment
disp('Processing extracellular compartment...')
[uwt_ex, dt_ex_total, dt_ex_all] = compute_compartment_tensors(V, V2, V3, fp, dt_ex_unit, wt_ex_unit);

% Combine compartments
dt_total = vf * dt_in_total + (1-vf) * dt_ex_total;

% Reshape 1D total diffusion tensor back to 3x3
DT = [dt_total(1), dt_total(4), dt_total(5);...
    dt_total(4), dt_total(2), dt_total(6);...
    dt_total(5), dt_total(6), dt_total(3)];

% Combine intrinsic microscopic kurtosis (weighted by volume and mean diffusivity squared)
KTi = vf * (mdi/md)^2 * uwt_in + (1-vf) * (mde/md)^2 * uwt_ex;

% Concatenate the orientation probability fractions for both compartments
fp_both = [fp * vf; fp * (1-vf)];

% Concatenate the 3x3 diffusion tensors for both compartments
num_dirs = length(V);
dt_all = zeros(3, 3, num_dirs * 2);
dt_all(:, :, 1:num_dirs) = dt_in_all;
dt_all(:, :, (num_dirs + 1):end) = dt_ex_all;

% Calculate covariance tensor (CT) and Variance Kurtosis (KTv)
CT = sim_zt(dt_all, DT, fp_both');
KTv = Zsymmetric(CT) / (md^2);

fff2 = 10;

% 3D geometries
np=100;
theta=linspace(0,2*pi,np+1);
phi=linspace(0,2*pi,np+1);
[theta, phi]=meshgrid(theta, phi);

scal = (length(parula):-1:1)/length(parula);
mycmap = [scal', zeros(length(parula), 2); 0, 0, 0; parula];


Evector = [0, 0, 1;...
    0, 1, 0;...
    1, 0, 0];

[dt_total_r, KTi] = fun_rotate_tensors_DKI(dt_total, KTi, Evector);
[dt_total, KTv] = fun_rotate_tensors_DKI(dt_total, KTv, Evector);

figure('color', [1 1 1])
subplot(1, 4, 1)
rd = DirectionalDiff_2D(dt_total, theta, phi);
[xx, yy, zz]=sph2cart(theta, phi, rd);
surf(xx, yy, zz, rd, 'EdgeColor', 'none')
title('Dapp'), xlabel('x'), ylabel('y'), zlabel('z')
view(20,20)
caxis([-2 2])
xlim([-2 2])
ylim([-2 2])
zlim([-2 2])
axis square

subplot(1, 4, 2)
r = DirectionalKurt_2D(dt_total, KTi+KTv, theta, phi);
[xx, yy, zz]=sph2cart(theta, phi, r);
surf(xx, yy, zz, r, 'EdgeColor', 'none')
title('W'), xlabel('x'), ylabel('y'), zlabel('z')
view(20,20)
caxis([-4 4])
xlim([-4 4])
ylim([-4 4])
zlim([-4 4])
axis square
colormap(mycmap)

subplot(1, 4, 4)
r = DirectionalKurt_2D(dt_total, KTi, theta, phi);
[xx, yy, zz]=sph2cart(theta, phi, r);
surf(xx, yy, zz, r, 'EdgeColor', 'none')
title('Wu'), xlabel('x'), ylabel('y'), zlabel('z')
view(20,20)
caxis([-1 1])
xlim([-1 1])
ylim([-1 1])
zlim([-1 1])
axis square
colormap(mycmap)


subplot(1, 4, 3)
r = DirectionalKurt_2D(dt_total, KTv, theta, phi);
[xx, yy, zz]=sph2cart(theta, phi, r);
surf(xx, yy, zz, r, 'EdgeColor', 'none')
title('Wv'), xlabel('x'), ylabel('y'), zlabel('z')
view(20,20)
caxis([-4 4])
xlim([-4 4])
ylim([-4 4])
zlim([-4 4])
axis square
colormap(mycmap)


%% CTI metrics
pars(1) = 1;
pars(2:7) = dt_total;
pars(8:22) = KTi+KTv;
pars(23:43) = CT;
[ktotal, kaniso, kiso, kintra, Convert, MDd] = fun_resolve_kurtosis_singlevoxel(pars);


[MD, AD, RD, FA, evec1, evec2, evec3, L2, L3] = dti_metrics(dt_total, 1);
[MKi, AKi, RKi] = dki_tensor_metrics(dt_total, KTi, 1);
[MKv, AKv, RKv] = dki_tensor_metrics(dt_total, KTv, 1);

fff =12;

figure('color', [1 1 1])
hAxes = subplot(1, 2, 1);
bar(1, ktotal, 'FaceColor',[30 30 30]/255)
hold on
bar(2, kaniso, 'FaceColor',[190 100 58]/255)
bar(3, kiso, 'FaceColor',[233 184 98]/255)
bar(4, kintra, 'FaceColor',[98 184 233]/255) %[68 154 217]/255
set(gca,'XTick',[1, 2, 3, 4]);
xticklabels({'$$K_{total}$$', '$$K_{aniso}$$','$$K_{iso}$$', '$${\mu}K$$'})
hAxes.TickLabelInterpreter = 'latex';
hAxes.FontSize = fff;
ylim([0, 3.5])



fig = figure('color', [1 1 1], 'Units', 'centimeters', ...
    'Position', [1 1 21 8]);

set(fig, 'PaperPositionMode', 'auto');
set(fig, 'PaperOrientation', 'portrait');
set(groot, 'DefaultTextFontName', 'Arial');


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subplot(1, 4, 1)
% this subplot will be used later for log norm distribution
axis off

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% CTI metrics
hAxes=subplot(1, 4, 2);

bar(1, ktotal, 'FaceColor',[30 30 30]/255)
hold on
bar(2, kaniso, 'FaceColor',[190 100 58]/255)
bar(3, kiso, 'FaceColor',[233 184 98]/255)
bar(4, kintra, 'FaceColor',[98 184 233]/255) %[68 154 217]/255
set(gca,'XTick',[0.8, 1.9, 3.1, 4.2]);

xticklabels({'$$\overline{K}_t$$', '$$\overline{K}_{ani}$$',...
    '$$\overline{K}_{iso}$$', '$$\overline{K}_{\mu}$$'})

xtickangle(hAxes, 0);
hAxes.TickLabelInterpreter = 'latex';
hAxes.FontSize = fff;
ylim([0, 2])
xlim([0.5, 4.5])
title({'B) Previous CTI'}, 'fontsize', fff2)
axis square

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Total Kurtosis metrics
subplot_tight(1, 6, 4)

r_kurt = DirectionalKurt_2D(dt_total, KTi+KTv, theta, phi);
[xx, yy, zz]=sph2cart(theta, phi, r_kurt);
surf(xx, yy, zz, r_kurt, 'EdgeColor', 'none')
view(30,30)
grid off
axis equal
caxis([-3.2 3.2])
xlim([-3.2 3.2])
ylim([-3.2 3.2])
zlim([-3.2 3.2])
title({'C) Total-Kurtosis';'(DKI)'}, 'fontsize', fff2)
text(1, 0, 3.2, '$$K_t(n)$$', 'Interpreter','latex', 'fontsize', fff2)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Variance Kurtosis metrics
subplot_tight(1, 6, 5)
r_kurt = DirectionalKurt_2D(dt_total, KTv, theta, phi);
[xx, yy, zz]=sph2cart(theta, phi, r_kurt);
surf(xx, yy, zz, r_kurt, 'EdgeColor', 'none')
view(30, 30)
grid off
axis equal
caxis([-3.2 3.2])
xlim([-3.2 3.2])
ylim([-3.2 3.2])
zlim([-3.2 3.2])

title({'D) Var-Kurtosis';'(new CTI)'}, 'fontsize', fff2)
text(1, 0, 3.2, '$$K_v(n)$$', 'Interpreter','latex', 'fontsize', fff2)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Micro Kurtosis metrics
subplot_tight(1, 6, 6)
r_kurt = DirectionalKurt_2D(dt_total, KTi, theta, phi);
[xx, yy, zz]=sph2cart(theta, phi, r_kurt);
surf(xx, yy, zz, r_kurt, 'EdgeColor', 'none')
view(30,30)
grid off
axis equal
caxis([-1 1])
xlim([-1 1])
ylim([-1 1])
zlim([-1 1])
title({'E) Micro-Kurtosis';'(new CTI)'}, 'fontsize', fff2)
text(1/3, 0, 3.2/3, '$$K_\mu(n)$$', 'Interpreter','latex', 'fontsize', fff2)


%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Log norm distribution
%% Running for two cases
%% case 1 - log-normal distribution parameters for Corpus Callosum Genu (m=3.0, std=0.6)
%% case 2 - log-normal distribution parameters for Spinal Cord (m=3.0, std=0.6)
fig = figure('color', [1 1 1], 'Units', 'centimeters', ...
    'Position', [1 1 21 11]);

set(fig, 'PaperPositionMode', 'auto');
set(fig, 'PaperOrientation', 'portrait');
set(groot, 'DefaultTextFontName', 'Arial');
for cases=1:2
    
    hAxes = subplot(2, 4, (cases-1)*4+1);
    
    if cases==1
        m = 1; % CCg
        v = 0.5^2;
    else
        m = 3; %% SC
        v = 0.6^2;
    end
    
    
    % Convert normal mean/var to log-normal mu/sigma
    mu = log((m^2)/sqrt(v+m^2));
    sigma = sqrt(log(v/(m^2)+1));
    
    % Define range of radii to simulate
    if cases==1
        r_sampled = 0.1:0.1:4;
    else
        r_sampled = 1.5:0.1:5.5;
    end
    
    % Calculate probability density (Count distribution)
    P_count = lognpdf(r_sampled, mu, sigma);
    
    % Volume weighting (signal fraction is proportional to cross-sectional area R^2)
    W_vol = P_count .* (r_sampled.^2);
    W_vol = W_vol / sum(W_vol); % Normalize weights to sum to 1
    
    bar(r_sampled, P_count, 1, 'FaceColor', [117 11 45]/255, 'EdgeColor', 'none');
    
    % Apply your specific formatting
    title(['A',num2str(cases),') Distribution radii'], 'fontsize', fff2, 'HorizontalAlignment', 'center')
    xlabel('$$r (\mu m)$$', 'Interpreter','latex', 'fontsize', fff)
    ylabel('$$P(r)$$', 'Interpreter','latex', 'fontsize', fff)
    axis square
    
    %%%%%%%%
    %%% Calculating distribution ensemble tensors
    %%%%%%%
    
    disp('Calculating distribution ensemble tensors...')
    
    num_radii = length(r_sampled);
    num_dirs = length(V);
    
    % 1. Pre-calculate the global Mean Diffusivity (md_global)
    mdi_r_array = zeros(1, num_radii);
    for ri = 1:num_radii
        [~, ~, mdi_r_array(ri)] = compute_cylinder_unit_tensors(10, r_sampled(ri), intraD/100, 20);
    end
    md_global = vf * sum(W_vol .* mdi_r_array) + (1-vf) * mde;
    
    % 2. Initialize arrays for the distributed compartments
    % Size: (orientations * radii) + (orientations for extracellular)
    total_compartments = num_dirs * num_radii + num_dirs;
    
    dt_all_dist = zeros(3, 3, total_compartments);
    fp_both_dist = zeros(total_compartments, 1);
    
    dt_in_total_dist = zeros(1, 6);
    KTi_dist = zeros(1, 15);
    
    idx = 1;
    
    % 3. Loop over the radii distribution
    for ri = 1:num_radii
        r_val = r_sampled(ri);
        
        % Get unit tensors for this specific radius
        [dt_unit_r, wt_unit_r, mdi_r] = compute_cylinder_unit_tensors(10, r_val, intraD/100, 20);
        
        % Apply the Bingham ODF orientation
        [uwt_in_r, dt_total_r, dt_all_r] = compute_compartment_tensors(V, V2, V3, fp, dt_unit_r, wt_unit_r);
        
        % Accumulate total intra diffusion (weighted by volume)
        dt_in_total_dist = dt_in_total_dist + W_vol(ri) * dt_total_r;
        
        % Accumulate microscopic kurtosis (weighted by volume fraction AND (MD_i/MD_global)^2 )
        KTi_dist = KTi_dist + (vf * W_vol(ri)) * (mdi_r / md_global)^2 * uwt_in_r;
        
        % Append the 3x3 matrices to the massive covariance array
        end_idx = idx + num_dirs - 1;
        dt_all_dist(:, :, idx:end_idx) = dt_all_r;
        
        % The probability weight for these specific tensors
        fp_both_dist(idx:end_idx) = fp * vf * W_vol(ri);
        
        idx = end_idx + 1;
    end
    
    % 4. Add the extracellular compartment to the ensemble
    KTi_dist = KTi_dist + (1-vf) * (mde / md_global)^2 * uwt_ex;
    dt_total_dist = vf * dt_in_total_dist + (1-vf) * dt_ex_total;
    
    end_idx = idx + num_dirs - 1;
    dt_all_dist(:, :, idx:end_idx) = dt_ex_all;
    fp_both_dist(idx:end_idx) = fp * (1-vf);
    
    % 5. Format total diffusion tensor and compute Variance Kurtosis
    DT_dist = [dt_total_dist(1), dt_total_dist(4), dt_total_dist(5);...
        dt_total_dist(4), dt_total_dist(2), dt_total_dist(6);...
        dt_total_dist(5), dt_total_dist(6), dt_total_dist(3)];
    
    CT_dist = sim_zt(dt_all_dist, DT_dist, fp_both_dist');
    KTv_dist = Zsymmetric(CT_dist) / (md_global^2);
    
    
    %% 6. Plotting the Distributed Results
    
    % Resolve scalar metrics for the bar chart
    pars_dist(1) = 1;
    pars_dist(2:7) = dt_total_dist;
    pars_dist(8:22) = KTi_dist + KTv_dist;
    pars_dist(23:43) = CT_dist;
    [ktotal_d, kaniso_d, kiso_d, kintra_d, ~, ~] = fun_resolve_kurtosis_singlevoxel(pars_dist);
    
    % Rotate to align with Z-axis for 3D plotting
    [~, KTi_dist_r] = fun_rotate_tensors_DKI(dt_total_dist, KTi_dist, Evector);
    [dt_total_dist_r, KTv_dist_r] = fun_rotate_tensors_DKI(dt_total_dist, KTv_dist, Evector);
    
    % Plot B: CTI Bar Chart
    hAxes = subplot(2, 4, (cases-1)*4+2);
    bar(1, ktotal_d, 'FaceColor',[30 30 30]/255); hold on;
    bar(2, kaniso_d, 'FaceColor',[190 100 58]/255);
    bar(3, kiso_d, 'FaceColor',[233 184 98]/255);
    bar(4, kintra_d, 'FaceColor',[98 184 233]/255);
    set(gca,'XTick',[0.8, 1.9, 3.1, 4.2]);
    xticklabels({'$$\overline{K}_t$$', '$$\overline{K}_{ani}$$', '$$\overline{K}_{iso}$$', '$$\overline{K}_{\mu}$$'})
    xtickangle(hAxes, 0);
    hAxes.TickLabelInterpreter = 'latex';
    hAxes.FontSize = fff;
    ylim([0, 2]); xlim([0.5, 4.5]);
    title(['B',num2str(cases),') Distributed CTI'], 'fontsize', fff2)
    axis square
    
    % Plot C: Total Kurtosis
    subplot_tight(2, 6, (cases-1)*6+4)
    r_kurt = DirectionalKurt_2D(dt_total_dist_r, KTi_dist_r + KTv_dist_r, theta, phi);
    [xx, yy, zz] = sph2cart(theta, phi, r_kurt);
    surf(xx, yy, zz, r_kurt, 'EdgeColor', 'none'); view(30,30); grid off; axis equal;
    caxis([-3.2 3.2]); xlim([-3.2 3.2]); ylim([-3.2 3.2]); zlim([-3.2 3.2]);
    title({['C',num2str(cases),') Total-Kurtosis'];'(DKI)'}, 'fontsize', fff2)
    text(1, 0, 3.2, '$$K_t(n)$$', 'Interpreter','latex', 'fontsize', fff2)
    
    % Plot D: Variance Kurtosis
    subplot_tight(2, 6, (cases-1)*6+5)
    rt = DirectionalKurt_2D(dt_total_dist_r, KTv_dist_r, theta, phi);
    [xx, yy, zz] = sph2cart(theta, phi, rt);
    surf(xx, yy, zz, rt, 'EdgeColor', 'none'); view(30, 30); grid off; axis equal;
    caxis([-3.2 3.2]); xlim([-3.2 3.2]); ylim([-3.2 3.2]); zlim([-3.2 3.2]);
    title({['D',num2str(cases),') Var-Kurtosis'];'(new CTI)'}, 'fontsize', fff2)
    text(1, 0, 3.2, '$$K_v(n)$$', 'Interpreter','latex', 'fontsize', fff2)
    
    % Plot E: Micro Kurtosis
    subplot_tight(2, 6, (cases-1)*6+6)
    r_mu = DirectionalKurt_2D(dt_total_dist_r, KTi_dist_r, theta, phi);
    [xx, yy, zz] = sph2cart(theta, phi, r_mu);
    surf(xx, yy, zz, r_mu, 'EdgeColor', 'none'); view(30,30); grid off; axis equal;
    caxis([-1 1]); xlim([-1 1]); ylim([-1 1]); zlim([-1 1]);
    title({['E',num2str(cases),') Micro-Kurtosis'];'(new CTI)'}, 'fontsize', fff2)
    text(1/3, 0, 3.2/3, '$$K_\mu(n)$$', 'Interpreter','latex', 'fontsize', fff2)
    
    disp('Done!')
end

print(fig, 'figureS5.pdf', '-dpdf', '-painters');
%print(fig, 'figureS5.eps', '-depsc', '-painters');
fprintf('Figures saved successfully!\n');

function [uwt_comp, dt_total_comp, dt_all_comp] = compute_compartment_tensors(V, V2, V3, fp, dt_unit, wt_unit)
% Helper function to rotate and accumulate tensors for a specific compartment
num_dirs = size(V, 1);

uwt_comp = 0;
dt_total_comp = 0;
dt_all_comp = zeros(3, 3, num_dirs);

for vi = 1:num_dirs
    Xa = V(vi, :);
    Xb = V2(vi, :);
    Xc = V3(vi, :);
    
    Evector = [Xa', Xb', Xc'];
    
    % Rotate the base tensors according to the Bingham ODF orientation
    [dt_rep, wt_rep] = fun_rotate_tensors_DKI(dt_unit, wt_unit, Evector);
    
    % Accumulate weighted total tensors
    uwt_comp = uwt_comp + fp(vi) * wt_rep;
    dt_total_comp = dt_total_comp + fp(vi) * dt_rep;
    
    % Map 1D diffusion tensor format back into 3x3 matrices for covariance steps
    dt_all_comp(1, 1, vi) = dt_rep(1);
    dt_all_comp(1, 2, vi) = dt_rep(4);
    dt_all_comp(1, 3, vi) = dt_rep(5);
    
    dt_all_comp(2, 1, vi) = dt_rep(4);
    dt_all_comp(2, 2, vi) = dt_rep(2);
    dt_all_comp(2, 3, vi) = dt_rep(6);
    
    dt_all_comp(3, 1, vi) = dt_rep(5);
    dt_all_comp(3, 2, vi) = dt_rep(6);
    dt_all_comp(3, 3, vi) = dt_rep(3);
end
end

function [dt_unit, wt_unit, mdi] = compute_cylinder_unit_tensors(diff_time, radius, d0, N_roots)
% Calculates the 1D diffusion and kurtosis tensors for a perfectly
% straight, impermeable cylinder of a given radius.

% Get radial metrics from the exact analytical Bessel solution
[RDi, RKi, ~] = diffusion_in_cylinder(diff_time, radius, N_roots, d0);

% Axial diffusivity is unrestricted (matches intrinsic diffusivity d0)
ADi = d0;

% Mean Diffusivity
mdi = (ADi + 2 * RDi) / 3;

% Axial Kurtosis is 0 for unrestricted 1D Gaussian diffusion
AKi = 0;

% Mean Kurtosis
MKi = (ADi^2 * AKi + (8/3) * RDi^2 * RKi) / (5 * mdi^2);

% Construct the 1D Diffusion Tensor [Dxx, Dyy, Dzz, Dxy, Dxz, Dyz]
% (Note: Dimension 1 is set as the axial dimension)
dt_unit = [ADi, RDi, RDi, 0, 0, 0];

% Construct the 1D Kurtosis Tensor (15 elements)
wt_unit = zeros(1, 15);
wt_unit(1) = AKi * ADi^2 / (mdi^2);
wt_unit(2) = RKi * RDi^2 / (mdi^2);
wt_unit(3) = wt_unit(2);

wt_unit(8) = 5.0 / 4 * MKi - wt_unit(1) / 4 - 2.0 / 3 * wt_unit(2);
wt_unit(9) = wt_unit(8);
wt_unit(10) = wt_unit(2) / 3;
end