%close all
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



%% Full CTI

ODI = 0.12;

coc = 1/tan(pi/2*ODI);
angle = asin(sqrt(1/(2*coc)))/pi*180;
load('dirs1024.mat')

fp = bingham(V, 1, coc, coc, [0, 1, 0]', [0, 0, 1]');
fp = fp/sum(fp);

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

dt_total = DT;

np=100;
theta=linspace(0,2*pi,np+1);
phi=linspace(0,2*pi,np+1);
[theta, phi]=meshgrid(theta, phi);

scal = (length(parula):-1:1)/length(parula);
mycmap = [scal', zeros(length(parula), 2); 0, 0, 0; parula];


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
%save('3D_ktotal_bead', 'xx', 'yy', 'zz', 'r')

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
%save('3D_kmicro_bead', 'xx', 'yy', 'zz', 'r')



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


save('case3', 'dt_total', 'KTi', 'KTv', ...
    'ktotal', 'kaniso', 'kiso', 'kintra', ...
    'MKi','AKi', 'RKi', 'MKv','AKv', 'RKv')