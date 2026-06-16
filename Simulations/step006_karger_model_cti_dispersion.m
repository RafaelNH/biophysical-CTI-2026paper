close all
clear all
clc

fs = filesep;
addpath(['..', fs, 'DWIu_toolbox_v1.11', fs, 'DKI_toolbox/'])
addpath(['..', fs, 'DWIu_toolbox_v1.11', fs, 'DDE_toolbox/'])
addpath(['..', fs, 'DWIu_toolbox_v1.11', fs, 'Spherical_functions_and_directions/'])
load('dirs1024.mat')

pth = ['MC_results_from_step001', fs];

extraD = 150;
intraD = 250;

seed = 123;

Delta = 10; %ms
Tmix = 10; %ms

tex = 3:0.1:1000; % ms
k = 1./tex; % ms^(-1)

As = 0;
f = 0.7;
as = As*100;

load([pth,'tensors_as',num2str(as),'_intraD',num2str(intraD),'_seed',(num2str(seed))]);
load([pth,'tensors_as',num2str(as),'_f',num2str(f*100),'_extraD',num2str(extraD),'_seed',num2str(seed)], 'wt_ex', 'dt_ex');
        

pth = ['sims', fs];

AS = 0:0.1:0.6;
vf = 0.7;
extraD = 150;
intraD = 250;

SEEDS = 123;

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

AK = AK0 * hsde;
RK = RK0 * hsde;

AKv = AK0 * hdde;
RKv = RK0 * hdde;

AKm = AK0 * (hsde-hdde);
RKm = RK0 * (hsde-hdde);

figure('color', [1 1 1])
subplot(2, 3, 1)
plot(k, AK)
hold on
yline(AK0)
title('AK')
%ylim([0, 0.01])

subplot(2, 3, 4)
plot(k, RK)
hold on
yline(RK0)
title('RK')

subplot(2, 3, 2)
plot(k, AKv)
hold on
yline(AK0)
title('AKv')

subplot(2, 3, 5)
plot(k, RKv)
hold on
yline(RK0)
title('RKv')

subplot(2, 3, 3)
plot(k, AKm)
hold on
yline(AK0)
title('AKm')

subplot(2, 3, 6)
plot(k, RKm)
hold on
yline(RK0)
title('RKm')

figure('color', [1 1 1])
subplot(1,2,1)
plot(k, AKv./AK)
hold on 
plot(k, AKm./AK)

subplot(1,2,2)
plot(k, RKv./RK)
hold on 
plot(k, RKm./RK)

%% Full CTI

ODI = 0.12;

coc = 1/tan(pi/2*ODI);
load('dirs1024.mat')

fp = bingham(V, 1, coc, coc, [0, 1, 0]', [0, 0, 1]');
fp = fp/sum(fp);

angle = asin(sqrt(1/(2*coc)))/pi*180;

[ZT0, KT0, DT] = Zgenerator_directions(ADi, RDi, ADe, RDe, f, fp', V);
VF =  2 * (k*Delta - (1 - exp(-k*Delta))) ./ ((k*Delta).^2);
ZF = exp(-k*Tmix) .* ((1 - exp(-k*Delta)).^2) ./ ((k*Delta).^2);
S0 = 1;

pars(1) = S0;
pars(2:7) = DT;
pars(8:22) = KT0;
pars(23:43) = ZT0;

[DT, KT, KTi, KTv] = fun_resolve_kurtosis_tensors_svoxel(pars);

DTensor =  [DT(1) DT(4) DT(5); DT(4) DT(2) DT(6); DT(5) DT(6) DT(3)];
[MDref, FAref, l1, l2, l3, v1, v2, v3] = DTImetrics(DTensor);
ADref = l1;
RDref = (l2+l3)/2;
[ktotalref, kanisoref, kisoref, kintraref, Convert, MDd] = fun_resolve_kurtosis_singlevoxel(pars);
[MKref, AKref, RKref] = dki_tensor_metrics(DT, KT, 1);
[MKiref, AKiref, RKiref] = dki_tensor_metrics(DT, KTi, 1);
[MKvref, AKvref, RKvref] = dki_tensor_metrics(DT, KTv, 1);


for si=1:length(k)
    KT = KT0 * VF(si);
    ZT = ZT0 * ZF(si);
    pars(1) = S0;
    pars(2:7) = DT;
    pars(8:22) = KT;
    pars(23:43) = ZT;
    
        DTensor =  [DT(1) DT(4) DT(5); DT(4) DT(2) DT(6); DT(5) DT(6) DT(3)];
    [MD(si), FA, l1, l2, l3, v1, v2, v3] = DTImetrics(DTensor);
    AD(si) = l1;
    RD(si) = (l2+l3)/2;
    
    [DT, KT, KTi, KTv] = fun_resolve_kurtosis_tensors_svoxel(pars);
    
    [ktotal(si), kaniso(si), kiso(si), kintra(si), Convert, MDd] = fun_resolve_kurtosis_singlevoxel(pars);
    [MK(si), AK(si), RK(si)] = dki_tensor_metrics(DT, KT, 1);
    [MKi(si), AKi(si), RKi(si)] = dki_tensor_metrics(DT, KTi, 1);
    [MKv(si), AKv(si), RKv(si)] = dki_tensor_metrics(DT, KTv, 1);
end

figure('color', [1 1 1])
subplot(2, 3, 1)
plot(k, AK)
hold on
yline(AKref)
title('AK')
%ylim([0, 0.01])

subplot(2, 3, 4)
plot(k, RK)
hold on
yline(RKref)
title('RK')

subplot(2, 3, 2)
plot(k, AKv)
hold on
yline(AKref)
title('AKv')

subplot(2, 3, 5)
plot(k, RKv)
hold on
yline(RKref)
title('RKv')

subplot(2, 3, 3)
plot(k, AKi)
hold on
yline(AKiref)
title('AKm')

subplot(2, 3, 6)
plot(k, RKi)
hold on
yline(RKiref)
title('RKm')


figure('color', [1 1 1])
subplot(2, 3, 1)
plot(k, ktotal)
hold on
yline(ktotalref)
title('ktotal')
%ylim([0, 0.01])

subplot(2, 3, 4)
plot(k, kaniso)
hold on
yline(kanisoref)
title('kaniso')

subplot(2, 3, 2)
plot(k, kiso)
hold on
yline(kisoref)
title('kiso')

subplot(2, 3, 5)
plot(k, kintra)
hold on
yline(kintraref)
title('kintra')



figure('color', [1 1 1], 'position', [-0.3333 455.0000 971.3333 250.6667])
subplot(1,3,3)
plot(k, AKv./AK, 'g')
hold on 
plot(k, AKi./AK, 'b')
xlim([0 0.33])
legend('$$vW_\parallel/W_\parallel$$', '$${\mu}W_\parallel/W_\parallel$$', 'Interpreter', 'latex','fontsize', 10)
xlabel('$$r (/ms)$$', 'Interpreter','latex', 'fontsize', 14)

subplot(1,3,2)
plot(k, RKv./RK, 'g')
hold on 
plot(k, RKi./RK, 'b')
xlim([0 0.33])
legend('$$vW_\bot/W_\bot$$', '$${\mu}W_\bot/W_\bot$$', 'Interpreter', 'latex','fontsize', 10)
xlabel('$$r (/ms)$$', 'Interpreter','latex', 'fontsize', 14)

subplot(1,3,1)
plot(k, MKv./MK, 'g')
hold on 
plot(k, MKi./MK, 'b')
xlim([0 0.33])
legend('$$v\overline{W}/\overline{W}$$', '$${\mu}\overline{W}/\overline{W}$$', 'Interpreter', 'latex','fontsize', 10)
xlabel('$$r (/ms)$$', 'Interpreter','latex', 'fontsize', 14)

figure('color', [1 1 1], 'position', [-0.3333 455.0000 971.3333 250.6667])
subplot(1, 5, 2)
plot(k, ktotal, 'b')
hold on
plot(k, kaniso, 'g')
plot(k, kiso, 'r')
plot(k, kintra, 'cyan')
xlim([0 0.33])
legend('$$Ktotal$$', '$$Kaniso$$','$$Kiso$$','$${\mu}K$$', 'Location','northeast','Interpreter','latex', 'fontsize', 10)
ylim([0 4])
xlabel('$$r (/ms)$$', 'Interpreter','latex', 'fontsize', 10)

subplot(1, 5, 1)
plot(k, MD, 'b')
hold on
plot(k, RD, 'g')
plot(k, AD, 'r')
xlim([0 0.33])
ylim([0 2])
legend('$$\overline{D}$$', '$$D_\bot$$', '$$D_\parallel$$', 'Location','northeast','Interpreter','latex', 'fontsize', 10)
xlabel('$$r (/ms)$$', 'Interpreter','latex', 'fontsize', 10)

subplot(1, 5, 3)
plot(k, MK, 'b')
hold on
plot(k, RK, 'g')
plot(k, AK, 'r')
ylim([0 4])
xlim([0 0.33])
%legend('MK*MD^2', 'RK*RD^2', 'AK*AD^2', 'Location','northwest')
legend('$$\overline{W}$$', '$$W_\bot$$', '$$W_\parallel$$', 'Location','northeast','Interpreter','latex', 'fontsize', 10)
xlabel('$$r (/ms)$$', 'Interpreter','latex', 'fontsize', 10)

subplot(1, 5, 5)
plot(k, MKi, 'b')
hold on
plot(k, RKi, 'g')
plot(k, AKi, 'r')
ylim([0 4])
xlim([0 0.33])
%legend('uMK*MD^2', 'uRK*RD^2', 'uAK*AD^2', 'Location','northwest')
legend('$${\mu}\overline{W}$$', '$${\mu}W_\bot$$', '$${\mu}W_\parallel$$', 'Location','northeast','Interpreter','latex', 'fontsize', 10)
xlabel('$$r (/ms)$$', 'Interpreter','latex', 'fontsize', 10)

subplot(1, 5, 4)
plot(k, MKv, 'b')
hold on
plot(k, RKv, 'g')
ylim([0 4])
plot(k, AKv, 'r')
xlim([0 0.33])
%legend('vMK*MD^2', 'vRK*RD^2', 'vAK*AD^2', 'Location','northwest')
legend('$$v\overline{W}$$', '$$vW_\bot$$', '$$vW_\parallel$$', 'Location','northeast','Interpreter','latex', 'fontsize', 10)
xlabel('$$r (/ms)$$', 'Interpreter','latex', 'fontsize', 10)





figure('color', [1 1 1], 'position', [-0.3333 455.0000 971.3333 250.6667])
subplot(1,3,3)
plot(k, AKv./AK, 'g')
hold on 
plot(k, AKi./AK, 'b')
xlim([0 0.33])
legend('$$vW_\parallel/W_\parallel$$', '$${\mu}W_\parallel/W_\parallel$$', 'Interpreter', 'latex','fontsize', 10)
xlabel('$$r (/ms)$$', 'Interpreter','latex', 'fontsize', 14)

subplot(1,3,2)
plot(k, RKv./RK, 'g')
hold on 
plot(k, RKi./RK, 'b')
xlim([0 0.33])
legend('$$vW_\bot/W_\bot$$', '$${\mu}W_\bot/W_\bot$$', 'Interpreter', 'latex','fontsize', 10)
xlabel('$$r (/ms)$$', 'Interpreter','latex', 'fontsize', 14)

subplot(1,3,1)
plot(k, MKv./MK, 'g')
hold on 
plot(k, MKi./MK, 'b')
xlim([0 0.33])
legend('$$v\overline{W}/\overline{W}$$', '$${\mu}\overline{W}/\overline{W}$$', 'Interpreter', 'latex','fontsize', 10)
xlabel('$$r (/ms)$$', 'Interpreter','latex', 'fontsize', 14)

figure('color', [1 1 1], 'position', [-0.3333 455.0000 971.3333 250.6667])

subplot(1, 3, 1)
plot(k, RK)
hold on
plot(k, AK)
ylim([0 4])
xlim([0 0.33])
%legend('MK*MD^2', 'RK*RD^2', 'AK*AD^2', 'Location','northwest')
legend('$$W_\bot$$', '$$W_\parallel$$', 'Location','northeast','Interpreter','latex', 'fontsize', 10)
xlabel('$$r (/ms)$$', 'Interpreter','latex', 'fontsize', 10)

subplot(1, 3, 3)
plot(k, RKi)
hold on
plot(k, AKi)
ylim([0 4])
xlim([0 0.33])
%legend('uMK*MD^2', 'uRK*RD^2', 'uAK*AD^2', 'Location','northwest')
legend('$${\mu}W_\bot$$', '$${\mu}W_\parallel$$', 'Location','northeast','Interpreter','latex', 'fontsize', 10)
xlabel('$$r (/ms)$$', 'Interpreter','latex', 'fontsize', 10)

subplot(1, 3, 2)
plot(k, RKv)
hold on
ylim([0 4])
plot(k, AKv)
xlim([0 0.33])
%legend('vMK*MD^2', 'vRK*RD^2', 'vAK*AD^2', 'Location','northwest')
legend('$$vW_\bot$$', '$$vW_\parallel$$', 'Location','northeast','Interpreter','latex', 'fontsize', 10)
xlabel('$$r (/ms)$$', 'Interpreter','latex', 'fontsize', 10)

r = k;
save('sim_kurt_vs_r', 'r',  'RKi', 'AKi', 'RKv', 'AKv', 'RK', 'AK')
