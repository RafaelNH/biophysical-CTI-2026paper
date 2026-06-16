close all
clear all
clc

fs = filesep;
addpath(['..', fs, 'DWIu_toolbox_v1.11', fs, 'DKI_toolbox/'])
addpath(['..', fs, 'DWIu_toolbox_v1.11', fs, 'DDE_toolbox/'])
addpath(['..', fs, 'DWIu_toolbox_v1.11', fs, 'Spherical_functions_and_directions/'])

%cases = 1; % intra/extra no beading
%cases = 2; % intra/extra with beading

for cases = 1:2

if cases == 1
    As = 0;
    vf = 0.7;
else
    As = 0.6;
    vf = 0.7;
end

extraD = 150;
intraD = 250;
as = As*100;
seed = 123;

pth = ['MC_results_from_step001', fs];
load([pth,'tensors_as',num2str(as),'_intraD',num2str(intraD),'_seed',(num2str(seed))]);
load([pth,'tensors_as',num2str(as),'_f',num2str(vf*100),'_extraD',num2str(extraD),'_seed',(num2str(seed))]);

% Change x axis to y axis of a single replica
Evector = [0, 0, 1;...
    0, 1, 0;...
    1, 0, 0];

mdi = mean(dt_in(1:3));
mde = mean(dt_ex(1:3));
md = vf * mdi + (1-vf) * mde;

kiso_check = (vf*(mdi - md)^2 + (1-vf)*(mde-md)^2) / (md^2) * 3;

[dt_in_unit, wt_in_unit] = fun_rotate_tensors_DKI(dt_in, wt_in, Evector);
[dt_ex_unit, wt_ex_unit] = fun_rotate_tensors_DKI(dt_ex, wt_ex, Evector);


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

nbins = 5;
fp = fp/(sum(fp)*nbins);
dt_in_all = zeros(3, 3, length(V)*nbins);
dt_ex_all = zeros(3, 3, length(V)*nbins);
fp_all = zeros(length(V)*nbins, 1);

uwt_in = 0;
uwt_ex = 0;
dt_in_total = 0;
dt_ex_total = 0;
ai = 0;
for bi = 1:nbins
    for vi = 1:size(V, 1)
        Xa = V(vi, :);
        Xb = V2(vi, :);
        Xc = V3(vi, :);
        
        Evector=[Xa',Xb',Xc'];
        
        [dt_in_rep, wt_in_rep] = fun_rotate_tensors_DKI(dt_in_unit, wt_in_unit, Evector);
        [dt_ex_rep, wt_ex_rep] = fun_rotate_tensors_DKI(dt_ex_unit, wt_ex_unit, Evector);
        
        uwt_in = uwt_in + fp(vi) * wt_in_rep; % since md is only different between intra and extra we can ignore the md correction factor
        uwt_ex = uwt_ex + fp(vi) * wt_ex_rep;
        dt_in_total = dt_in_total + fp(vi) * dt_in_rep;
        dt_ex_total = dt_ex_total + fp(vi) * dt_ex_rep;
        
        ai = ai + 1;
        dt_in_all(1, 1, ai) = dt_in_rep(1);
        dt_in_all(1, 2, ai) = dt_in_rep(4);
        dt_in_all(1, 3, ai) = dt_in_rep(5);
        dt_in_all(2, 1, ai) = dt_in_rep(4);
        dt_in_all(2, 2, ai) = dt_in_rep(2);
        dt_in_all(2, 3, ai) = dt_in_rep(6);
        dt_in_all(3, 1, ai) = dt_in_rep(5);
        dt_in_all(3, 2, ai) = dt_in_rep(6);
        dt_in_all(3, 3, ai) = dt_in_rep(3);
        
        dt_ex_all(1, 1, ai) = dt_ex_rep(1);
        dt_ex_all(1, 2, ai) = dt_ex_rep(4);
        dt_ex_all(1, 3, ai) = dt_ex_rep(5);
        dt_ex_all(2, 1, ai) = dt_ex_rep(4);
        dt_ex_all(2, 2, ai) = dt_ex_rep(2);
        dt_ex_all(2, 3, ai) = dt_ex_rep(6);
        dt_ex_all(3, 1, ai) = dt_ex_rep(5);
        dt_ex_all(3, 2, ai) = dt_ex_rep(6);
        dt_ex_all(3, 3, ai) = dt_ex_rep(3);
        fp_all(ai) = fp(vi);
    end
    disp(bi)
end

dt_total = vf * dt_in_total + (1-vf) * dt_ex_total;

DT = [dt_total(1), dt_total(4), dt_total(5);...
    dt_total(4), dt_total(2), dt_total(6);...
    dt_total(5), dt_total(6), dt_total(3)];

KTi = vf * (mdi/md)^2 * uwt_in + (1-vf) * (mde/md)^2 * uwt_ex;

fp_both = [fp_all*vf; fp_all*(1-vf)];

dt_all = zeros(3, 3, length(V)*nbins*2);
dt_all(:, :, 1:(length(V)*nbins)) = dt_in_all;
dt_all(:, :, (1+length(V)*nbins):(length(V)*nbins*2)) = dt_ex_all;

CT = sim_zt(dt_all, DT, fp_both');

KTv = Zsymmetric(CT)/(md^2);


% Test DP-DKI metrics
% KT = KTi+KTv;
% Wbar = (1/5)*(KT(1) + KT(2) + KT(3) + 2*KT(10) + 2*KT(11) + 2*KT(12));
% Wtbar = (1/8)*(KT(1) + KT(2) + KT(3) + 2*KT(10) + 2*KT(11) + 2*KT(12)+...
%                CT(1)/(md^2) + CT(2)/(md^2) + CT(3)/(md^2) + 2*CT(10)/(md^2) + 2*CT(11)/(md^2) + 2*CT(12)/(md^2));
%      
% deltaW = Wbar - Wtbar;
% Wplus = 1/10 * (KT(1)   +   KT(2) + KT(3) +   2*KT(10) + 2*KT(11) + 2*KT(12)+...
%                 3*CT(1)/(md^2) + 3*CT(2)/(md^2) + 3*CT(3)/(md^2) + 2*CT(10)/(md^2) + 2*CT(11)/(md^2) + 2*CT(12)/(md^2)+...
%                 4*CT(16)/(md^2) + 4*CT(17)/(md^2) + 4*CT(18)/(md^2));
% Wminus = Wplus; % Because I am considering that S tensor vanishes.
% dWex = Wbar - 0.5*(Wplus+Wminus)

% 3D geometries
np=100;
theta=linspace(0,2*pi,np+1);
phi=linspace(0,2*pi,np+1);
[theta, phi]=meshgrid(theta, phi);

scal = (length(parula):-1:1)/length(parula);
mycmap = [scal', zeros(length(parula), 2); 0, 0, 0; parula];

% intra and extra compartments
% Change x axis to y axis
Evector = [0, 0, 1;...
    0, 1, 0;...
    1, 0, 0];

[dt_in_total, uwt_in] = fun_rotate_tensors_DKI(dt_in_total, uwt_in, Evector);
[dt_ex_total, uwt_ex] = fun_rotate_tensors_DKI(dt_ex_total, uwt_ex, Evector);


figure('color', [1 1 1])
subplot(1, 4, 1)
rd = DirectionalDiff_2D(dt_in_total, theta, phi);
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
r = DirectionalKurt_2D(dt_in_total, uwt_in, theta, phi);
[xx, yy, zz]=sph2cart(theta, phi, r);
surf(xx, yy, zz, r, 'EdgeColor', 'none')
title('KT'), xlabel('x'), ylabel('y'), zlabel('z')
view(20,20)
caxis([-0.5 0.5])
xlim([-0.5 0.5])
ylim([-0.5 0.5])
zlim([-0.5 0.5])
axis square
colormap(mycmap)

subplot(1, 4, 3)
rd = DirectionalDiff_2D(dt_ex_total, theta, phi);
[xx, yy, zz]=sph2cart(theta, phi, rd);
surf(xx, yy, zz, rd, 'EdgeColor', 'none')
title('Dapp'), xlabel('x'), ylabel('y'), zlabel('z')
view(20,20)
caxis([-2 2])
xlim([-2 2])
ylim([-2 2])
zlim([-2 2])
axis square

subplot(1, 4, 4)
r = DirectionalKurt_2D(dt_ex_total, uwt_ex, theta, phi);
[xx, yy, zz]=sph2cart(theta, phi, r);
surf(xx, yy, zz, r, 'EdgeColor', 'none')
title('KT'), xlabel('x'), ylabel('y'), zlabel('z')
view(20,20)
caxis([-0.5 0.5])
xlim([-0.5 0.5])
ylim([-0.5 0.5])
zlim([-0.5 0.5])
axis square
colormap(mycmap)


% combined compartments
% Change x axis to y axis
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
%save('3D_kvar_bead', 'xx', 'yy', 'zz', 'r')


%% CTI metrics
pars(1) = 1;
pars(2:7) = dt_total;
pars(8:22) = KTi+KTv;
pars(23:43) = CT;
[ktotal, kaniso, kiso, kintra, Convert, MDd] = fun_resolve_kurtosis_singlevoxel(pars);


[MD, AD, RD, FA, evec1, evec2, evec3, L2, L3] = dti_metrics(dt_total, 1);
[MKi, AKi, RKi] = dki_tensor_metrics(dt_total, KTi, 1);
[MKv, AKv, RKv] = dki_tensor_metrics(dt_total, KTv, 1);

save(['case', num2str(cases)], 'dt_total', 'KTi', 'KTv', ...
    'ktotal', 'kaniso', 'kiso', 'kintra', ...
    'AKi', 'RKi', 'AKv', 'RKv')

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

hAxes = subplot(1, 2, 2);
bar(1, RKv, 'FaceColor',[117 11 45]/255) %[117 11 45]/255
hold on
bar(2, AKv, 'FaceColor',[217 154 68]/255) %[233 184 98]/255
bar(3, RKi, 'FaceColor',[45 11 117]/255) %[45 11 117]/255
bar(4, AKi, 'FaceColor',[68 154 217]/255) %[98 184 233]/255
set(gca,'XTick',[1, 2, 3, 4]);
xticklabels({'$$vW_\bot$$', '$$vW_\parallel$$','$${\mu}W_\bot$$', '$${\mu}W_\parallel$$'})
hAxes.TickLabelInterpreter = 'latex';
hAxes.FontSize = fff;
ylim([0, 3.5])

save(['case', num2str(cases)], 'dt_total', 'KTi', 'KTv', ...
    'ktotal', 'kaniso', 'kiso', 'kintra', ...
    'MKi','AKi', 'RKi', 'MKv','AKv', 'RKv')

end