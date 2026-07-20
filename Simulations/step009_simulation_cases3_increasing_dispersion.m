close all
clear all
clc

fs = filesep;
addpath(['..', fs, 'DWIu_toolbox_v1.11', fs, 'DKI_toolbox/'])
addpath(['..', fs, 'DWIu_toolbox_v1.11', fs, 'DDE_toolbox/'])
addpath(['..', fs, 'DWIu_toolbox_v1.11', fs, 'Spherical_functions_and_directions/'])

Delta = 10; %ms
Tmix = 10; %ms

tex = 33; % ms
k = 1./tex; % ms^(-1)

VF =  2 * (k*Delta - (1 - exp(-k*Delta))) ./ ((k*Delta).^2);
ZF = exp(-k*Tmix) .* ((1 - exp(-k*Delta)).^2) ./ ((k*Delta).^2);

As = 0;
f = 0.7;
extraD = 150;
intraD = 250;
as = As*100;
seed = 123;


pth = ['MC_results_from_step001', fs];
load([pth,'tensors_as',num2str(as),'_intraD',num2str(intraD),'_seed',(num2str(seed))]);
load([pth,'tensors_as',num2str(as),'_f',num2str(f*100),'_extraD',num2str(extraD),'_seed',(num2str(seed))]);


ADi = dt_in(3);
ADe = dt_ex(3);

RDi = (dt_in(1)+dt_in(2))/2;
RDe = (dt_ex(1)+dt_ex(2))/2;

AD = f*ADi + (1-f)*ADe;
RD = f*RDi + (1-f)*RDe;

D = (AD + 2*RD)/3;
Di = (ADi + 2*RDi)/3;
De = (ADe + 2*RDe)/3;

AK0 = 3*(f*(ADi - AD)^2 + (1-f)*(ADe - AD)^2)/(AD^2);
RK0 = 3*(f*(RDi - RD)^2 + (1-f)*(RDe - RD)^2)/(RD^2);

t1 = k*Delta;
T = 2*Delta + Tmix;
t = k*T;

hsde = 2./t1 - 2./(t1.^2) + 2*exp(-t1)./(t1.^2);
hdde = (exp(-t+2*t1)+exp(-t)-2*exp(-t+t1))./(t1.^2);


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
    fp = fp/sum(fp);
    [ZT0, KT0, DT] = Zgenerator_directions(ADi, RDi, ADe, RDe, f, fp', V);
    S0 = 1;

    KT = KT0 * VF;
    ZT = ZT0 * ZF;
    pars(1) = S0;
    pars(2:7) = DT;
    pars(8:22) = KT;
    pars(23:43) = ZT;
    
    DTensor =  [DT(1) DT(4) DT(5); DT(4) DT(2) DT(6); DT(5) DT(6) DT(3)];
    [MD(op_ind), FA(op_ind), l1, l2, l3, v1, v2, v3] = DTImetrics(DTensor);
    AD(op_ind) = l1;
    RD(op_ind) = (l2+l3)/2;
    
    [DT, KT, KTi, KTv] = fun_resolve_kurtosis_tensors_svoxel(pars);
    
    [ktotal(op_ind), kaniso(op_ind), kiso(op_ind), kintra(op_ind), Convert, MDd] = fun_resolve_kurtosis_singlevoxel(pars);
    
    [MKi(op_ind), AKi(op_ind), RKi(op_ind)] = dki_tensor_metrics(DT, KTi, 1);
    [MKv(op_ind), AKv(op_ind), RKv(op_ind)] = dki_tensor_metrics(DT, KTv, 1);
    
    % Estimation based on MDE literature
    uFA2_cti = 15 * kaniso(op_ind) / (10 * kaniso(op_ind) + 12);
    OP_estimate_1(op_ind) = sqrt((3/(uFA2_cti)-2)/(3/(FA(op_ind)^2)-2));
    
    % I will use the estimation based on literature, however, I am testing
    % bellow an estimate that may be more accurate
    % Using this will be irrelevant as Kiso are typically low in brain data
    uFA2_cti = 15 * kaniso(op_ind) / (10 * kaniso(op_ind) + 4*kiso(op_ind) + 12);
    OP_estimate_2(op_ind) = sqrt((3/(uFA2_cti)-2)/(3/(FA(op_ind)^2)-2));
     
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

save('sim_kurt_vs_op_exchange')



