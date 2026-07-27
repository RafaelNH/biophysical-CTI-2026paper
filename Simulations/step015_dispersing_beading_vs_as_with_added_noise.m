close all
clear all
clc

fs = filesep;
addpath(['..', fs, 'DWIu_toolbox_v1.11', fs, 'DKI_toolbox/'])
addpath(['..', fs, 'DWIu_toolbox_v1.11', fs, 'DDE_toolbox/'])
addpath(['..', fs, 'DWIu_toolbox_v1.11', fs, 'Spherical_functions_and_directions/'])
load('dirs1024.mat')

% parameters for noise simulations
load('gtab')
SNR = 40;
nrep = 1000;
A = fun_CTI_desing_matrix(gtab.bval1(1:263), gtab.bval2(1:263), ...
    gtab.dir1(:, 1:263), gtab.dir2(:, 1:263));
piA = pinv(A);
sig = 1/SNR;
Nvol = size(A, 1);


pth = ['MC_results_from_step001', fs];

AS = 0:0.1:0.6;

% volume fraction
vf = 0.7;

extraD = 150;
intraD = 250;

SEEDS = 123;

ODI = 0.12;
k = 1/tan(pi/2*ODI);
angle = asin(sqrt(1/(2*k)))/pi*180;

% variables that allow to set intra or extra cellular kurtosis to zero
% (used for supplementary material S1), for main simulation please maintain
% these variables to false)
inside_kurtosis_only = false;
outside_kurtosis_only = false;

aiii = 0;
for As = AS
    aiii = aiii+1;
    as = As*100;
    siii = 0;
    for seed = SEEDS
        siii = siii + 1;
        %vf = 0.7;
        load([pth,'tensors_as',num2str(as),'_intraD',num2str(intraD),'_seed',(num2str(seed))]);
        load([pth,'tensors_as',num2str(as),'_f',num2str(vf*100),'_extraD',num2str(extraD),'_seed',num2str(seed)], 'wt_ex', 'dt_ex');
        % Change x axis to y axis of a single replica
        Evector = [0, 0, 1;...
            0, 1, 0;...
            1, 0, 0];
        
        if inside_kurtosis_only
            wt_ex = wt_ex*0;
        elseif outside_kurtosis_only
            wt_in = wt_in*0;
        end
        [dt_in_unit, wt_in_unit] = fun_rotate_tensors_DKI(dt_in, wt_in, Evector);
        [dt_ex_unit, wt_ex_unit] = fun_rotate_tensors_DKI(dt_ex, wt_ex, Evector);
        
        %vf = 0.7;
        
        mdi = mean(dt_in_unit(1:3));
        mde = mean(dt_ex_unit(1:3));
        md = vf * mdi + (1-vf) * mde;
        
        % Add dispersion
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
        
        %% CTI metrics
        pars(1) = 1;
        pars(2:7) = dt_total;
        pars(8:22) = KTi+KTv;
        pars(23:43) = CT;
        [ktotal(aiii, siii), kaniso(aiii, siii), kiso(aiii, siii), kintra(aiii, siii), Convert, MDd] = fun_resolve_kurtosis_singlevoxel(pars);
        
        [MD(aiii, siii), AD(aiii, siii), RD(aiii, siii), FA(aiii, siii), evec1, evec2, evec3, L2, L3] = dti_metrics(dt_total, 1);
        [MK(aiii, siii), AK(aiii, siii), RK(aiii, siii)] = dki_tensor_metrics(dt_total, KTi+KTv, 1);
        [MKi(aiii, siii), AKi(aiii, siii), RKi(aiii, siii)] = dki_tensor_metrics(dt_total, KTi, 1);
        [MKv(aiii, siii), AKv(aiii, siii), RKv(aiii, siii)] = dki_tensor_metrics(dt_total, KTv, 1);
        
        % Generate predicted signals
        X(43) = log(1);
        X(1:6) = dt_total;
        X(7:21) = (KTi+KTv) * MD(aiii, siii)^2;
        X(22:42) = CT;
        Spred = exp(A*X');
        
        % run for different noise instances
        for ri = nrep:-1:1
            GR = sig*randn(Nvol,1);
            GI = sig*randn(Nvol,1);
            S = sqrt((Spred/sqrt(2)+GR).^2+(Spred/sqrt(2)+GI).^2);
            
            X0=piA*log(S);
            
            pars(1) = exp(X0(43));
            pars(2:7) = X0(1:6);
            pars(8:22) = X0(7:21)/(mean(X0(1:3))^2);
            pars(23:43) = X0(22:42);
            
            [ktotal_noise(aiii, ri), kaniso_noise(aiii, ri), kiso_noise(aiii, ri), kintra_noise(aiii, ri), Convert, MDd] = fun_resolve_kurtosis_singlevoxel(pars);
        
            [DT_noise, KT_noise, KTi_noise, KTv_noise] = fun_resolve_kurtosis_tensors_svoxel(pars);
            
            [MD_noise(aiii, ri), AD_noise(aiii, ri), RD_noise(aiii, ri), FA_noise(aiii, ri), evec1, evec2, evec3, L2, L3] = dti_metrics(DT_noise, 1);
            [MK_noise(aiii, ri), AK_noise(aiii, ri), RK_noise(aiii, ri)] = dki_tensor_metrics(DT_noise, KTi_noise+KTv_noise, 1);
            [MKi_noise(aiii, ri), AKi_noise(aiii, ri), RKi_noise(aiii, ri)] = dki_tensor_metrics(DT_noise, KTi_noise, 1);
            [MKv_noise(aiii, ri), AKv_noise(aiii, ri), RKv_noise(aiii, ri)] = dki_tensor_metrics(DT_noise, KTv_noise, 1);
        
            
        end
         disp(['As ', num2str(as)])
        
    end
end

save('Noise_simulations')

