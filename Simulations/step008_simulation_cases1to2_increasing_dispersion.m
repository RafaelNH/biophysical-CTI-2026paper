close all
clear all
clc

fs = filesep;
addpath(['..', fs, 'DWIu_toolbox_v1.11', fs, 'DKI_toolbox/'])
addpath(['..', fs, 'DWIu_toolbox_v1.11', fs, 'DDE_toolbox/'])
addpath(['..', fs, 'DWIu_toolbox_v1.11', fs, 'Spherical_functions_and_directions/'])

%cases = 1; % intra/extra no beading
%cases = 2; % intra/extra with beading
cases = 2;

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

[mdi, adi, rdi, FAi, evec1, evec2, evec3, L2, L3] = dti_metrics(dt_in, 1);
[mde, ade, rde, FAe, evec1, evec2, evec3, L2, L3] = dti_metrics(dt_ex, 1);

VAReig = vf*((adi-mdi)^2 + 2*(rdi-mdi)^2)/3 + (1-vf)*((ade-mde)^2 + 2*(rde-mde)^2)/3;
kaniso_check = 6/5 * VAReig / (md^2);

uFA2_check = 3/2  * VAReig/(VAReig + md^2);
uFA2_cti = kaniso_check/(kaniso_check + 6/5) * 3/2;
[dt_in_unit, wt_in_unit] = fun_rotate_tensors_DKI(dt_in, wt_in, Evector);
[dt_ex_unit, wt_ex_unit] = fun_rotate_tensors_DKI(dt_ex, wt_ex, Evector);


load('dirs1024.mat')
disp('sampling bingham...')
OPs = [0.01:0.01:0.9];
kappas = OP2kappa(OPs);
angles = asin(sqrt(1./(2*kappas)))/pi*180;
ODIs = 2/pi * atan(1./kappas);

figure
subplot(2, 2, 1)
plot(OPs, kappas)
xlabel('OPs')
ylabel('k')

subplot(2, 2, 2)
plot(OPs, angles)
xlabel('OPs')
ylabel('FWHM')

subplot(2, 2, 3)
plot(OPs, ODIs)
xlabel('OPs')
ylabel('ODI')

V2 = V;
V3 = V;

op_ind = 0;
for k = kappas
    op_ind = op_ind + 1;
    
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
    
    nbins = 1;
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
    
    %% CTI metrics
    pars(1) = 1;
    pars(2:7) = dt_total;
    pars(8:22) = KTi+KTv;
    pars(23:43) = CT;
    [ktotal(op_ind), kaniso(op_ind), kiso(op_ind), kintra(op_ind), Convert, MDd] = fun_resolve_kurtosis_singlevoxel(pars);
    
    [MD(op_ind), AD(op_ind), RD(op_ind), FA(op_ind), evec1, evec2, evec3, L2, L3] = dti_metrics(dt_total, 1);
    
    
    % Estimation based on MDE literature
    uFA2_cti = 15 * kaniso(op_ind) / (10 * kaniso(op_ind) + 12);
    OP_estimate_1(op_ind) = sqrt((3/(uFA2_cti)-2)/(3/(FA(op_ind)^2)-2));
    
    % I will use the estimation based on literature, however, I am testing
    % bellow an estimate that may be more accurate
    % Using this will be irrelevant as Kiso are typically low in brain data
    uFA2_cti = 15 * kaniso(op_ind) / (10 * kaniso(op_ind) + 4*kiso(op_ind) + 12);
    OP_estimate_2(op_ind) = sqrt((3/(uFA2_cti)-2)/(3/(FA(op_ind)^2)-2));
   
    [MKi(op_ind), AKi(op_ind), RKi(op_ind)] = dki_tensor_metrics(dt_total, KTi, 1);
    [MKv(op_ind), AKv(op_ind), RKv(op_ind)] = dki_tensor_metrics(dt_total, KTv, 1);
    
    disp(['Done ', num2str(op_ind/length(OPs)*100), '%'])
end


fff =12;

% The following figure shows that measured OP from CTI parameters are close
% to OP ground truth values
figure('color', [1 1 1])
plot(OPs, OP_estimate_1, 'red')
hold on
plot(OPs, OP_estimate_2, 'green')
plot(OPs, OPs, '--')

% OP influence on kurtosis estimates

figure('color', [1 1 1])
hAxes = subplot(2, 2, 1);
plot(OPs, AKi, 'red')
hold on
plot(OPs, RKi, 'green')
legend('axial', 'radial')
xlabel('OP')
title('microscopic kurtosis')

hAxes = subplot(2, 2, 2);
plot(OPs, AKv, 'red')
hold on
plot(OPs, RKv, 'green')
legend('axial', 'radial')
title('variance kurtosis')
xlabel('OP')

hAxes = subplot(2, 2, 3);
plot(ODIs, AKi, 'red')
hold on
plot(ODIs, RKi, 'green')
legend('axial', 'radial')
xlabel('ODI')
title('microscopic kurtosis')

hAxes = subplot(2, 2, 4);
plot(ODIs, AKv, 'red')
hold on
plot(ODIs, RKv, 'green')
legend('axial', 'radial')
title('variance kurtosis')
xlabel('ODI')

ind = 66;

save('sim_kurt_vs_op_bead')



