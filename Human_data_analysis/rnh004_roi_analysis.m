close all
clear all
clc

fs = filesep;
addpath(['..', fs, 'NIFTI_toolbox'])
addpath(['..', fs, 'DWIu_toolbox_v1.11', fs, 'DDE_toolbox'])

Len = 10;

pth = 'C:\Users\rafae\Data\CTIhuman\sub-';
fold = 'CTImetrics_gs_c2';

roi_names = textread('JHU-labels.txt','%s');

% sort OP or RKu
fff = 14;
sort_M = 1;

%rois = [1:6, 7:2:47];
rois = [1, 3:5, 7, 15:2:37, 41, 45];
nvoxel_wm_mat = zeros(Len, length(rois));
for s = Len:-1:1
    if s < 10
        V = load_untouch_nii([pth, '0', num2str(s), fs, fold, fs, 'ROIs.nii']);
    else
        V = load_untouch_nii([pth, num2str(s), fs, fold, fs, 'ROIs.nii']);
    end
    WM=V.img;
    for ri = 1:length(rois)
        ri_ind = rois(ri);
        if ri_ind < 7
            nvoxel_wm_mat(s, ri) = sum(WM(:) == ri_ind);
        else
            nvoxel_wm_mat(s, ri) = sum(WM(:) == ri_ind) + sum(WM(:) == ri_ind+1);
        end
    end
end


% Prepare names
Nroi = 0; new_names = {};
for r = rois % unimodal ROIs (i.e. 6 first ROIs)
    Nroi = Nroi+1;
    if r < 7
        name = roi_names{r};
        name(find(name=='_')) = ' ';
    else
        name = roi_names{r}(1:(end-2));
    end
    name(find(name == '_')) = ' ';
    new_names{Nroi} = name;
end

% average number of voxels
mnv = mean(nvoxel_wm_mat, 1);
figure('color', [1 1 1]);
subplot(2,1,1),axis([0.5 27.5 0 1])
plot(mnv)
yline(100)
xlim([0.5, 27.5])
%set(gca, 'Xtick', 1:27, 'XTickLabels', {new_names})
subplot(2,1,2),axis([0.5 27.5 0 1])
for r = 1:Nroi
    t = text(r,0,new_names{r},'Rotation',90,'FontSize',10);
end
axis off

% template values
FA_wm_temp=zeros(1, Nroi);
uFA_wm_temp=zeros(1, Nroi);
MD_wm_temp=zeros(1, Nroi);
AD_wm_temp=zeros(1, Nroi);
RD_wm_temp=zeros(1, Nroi);
MK_wm_temp=zeros(1, Nroi);
AK_wm_temp=zeros(1, Nroi);
RK_wm_temp=zeros(1, Nroi);
MKv_wm_temp=zeros(1, Nroi);
AKv_wm_temp=zeros(1, Nroi);
RKv_wm_temp=zeros(1, Nroi);
MKu_wm_temp=zeros(1, Nroi);
AKu_wm_temp=zeros(1, Nroi);
RKu_wm_temp=zeros(1, Nroi);

load_templated

for ri = 1:length(rois)
    ri_ind = rois(ri);
    if ri_ind < 7
        WMroi = ROI_temp(:) == ri_ind;
    else
        WMroi = (ROI_temp(:) == ri_ind) + (ROI_temp(:) == ri_ind+1);
    end
    FA_wm_temp(ri) = mean(FA(WMroi>0));
    MD_wm_temp(ri) = mean(MD(WMroi>0));
    AD_wm_temp(ri) = mean(AD(WMroi>0));
    RD_wm_temp(ri) = mean(RD(WMroi>0));
    
    MK_wm_temp(ri) = mean(MK(WMroi>0));
    AK_wm_temp(ri) = mean(AK(WMroi>0));
    RK_wm_temp(ri) = mean(RK(WMroi>0));
    
    MKv_wm_temp(ri) = mean(MKv(WMroi>0));
    AKv_wm_temp(ri) = mean(AKv(WMroi>0));
    RKv_wm_temp(ri) = mean(RKv(WMroi>0));
    
    MKu_wm_temp(ri) = mean(MKi(WMroi>0));
    AKu_wm_temp(ri) = mean(AKi(WMroi>0));
    RKu_wm_temp(ri) = mean(RKi(WMroi>0));
    
    uFA_wm_temp(ri) = mean(uFA(WMroi>0));
    
end
c = 4;
x = 1/(3*c+6);
xs = (c+1)*x;
figure('color', [1 1 1])
subplot(2,1,1)
bar((1:Nroi)-xs, MD_wm_temp, c*x)
hold on
bar(1:Nroi, RD_wm_temp, c*x)
bar((1:Nroi)+xs, AD_wm_temp, c*x)
xlim([0.5, Nroi+0.5])
set(gca,'XTick',[],'YTick',[0 0.5 1 1.5 2]);
title('Diff')
subplot(4,1,3),axis([0.5 Nroi+0.5 0 1])
for r = 1:Nroi
    t = text(r,0,new_names{r},'Rotation',90,'FontSize',10);
end
axis off


figure('color', [1 1 1])
subplot(2,1,1)
bar((1:Nroi)-xs, MK_wm_temp, c*x)
hold on
bar(1:Nroi, RK_wm_temp, c*x)
bar((1:Nroi)+xs, AK_wm_temp, c*x)
xlim([0.5, Nroi+0.5])
set(gca,'XTick',[],'YTick',[0 0.5 1 1.5 2]);
title('Kurt')
subplot(4,1,3),axis([0.5 Nroi+0.5 0 1])
for r = 1:Nroi
    t = text(r,0,new_names{r},'Rotation',90,'FontSize',10);
end
axis off


figure('color', [1 1 1])
subplot(2,1,1)
bar((1:Nroi)-xs, MKv_wm_temp, c*x)
hold on
bar(1:Nroi, RKv_wm_temp, c*x)
bar((1:Nroi)+xs, AKv_wm_temp, c*x)
xlim([0.5, Nroi+0.5])
set(gca,'XTick',[],'YTick',[0 0.5 1 1.5 2]);
title('Var Kurt')
subplot(4,1,3),axis([0.5 Nroi+0.5 0 1])
for r = 1:Nroi
    t = text(r,0,new_names{r},'Rotation',90,'FontSize',10);
end
axis off

figure('color', [1 1 1])
subplot(2,1,1)
bar((1:Nroi)-xs, MKu_wm_temp, c*x)
hold on
bar(1:Nroi, RKu_wm_temp, c*x)
bar((1:Nroi)+xs, AKu_wm_temp, c*x)
xlim([0.5, Nroi+0.5])
set(gca,'XTick',[],'YTick',[0 0.5 1 1.5 2]);
title('Micro Kurt')
subplot(4,1,3),axis([0.5 Nroi+0.5 0 1])
for r = 1:Nroi
    t = text(r,0,new_names{r},'Rotation',90,'FontSize',10);
end
axis off


figure('color', [1 1 1])
subplot(2,1,1)
bar((1:Nroi)-xs, MKv_wm_temp./MK_wm_temp, c*x)
hold on
bar(1:Nroi, RKv_wm_temp./RK_wm_temp, c*x)
bar((1:Nroi)+xs, AKv_wm_temp./AK_wm_temp, c*x)
xlim([0.5, Nroi+0.5])
set(gca,'XTick',[],'YTick',[0 0.5 1]);
title('Kv/Kt')
subplot(4,1,3),axis([0.5 Nroi+0.5 0 1])
for r = 1:Nroi
    t = text(r,0,new_names{r},'Rotation',90,'FontSize',10);
end
axis off

figure('color', [1 1 1])
subplot(2,1,1)
bar((1:Nroi)-xs, FA_wm_temp, c*x)
hold on
bar(1:Nroi, uFA_wm_temp, c*x)
%bar((1:Nroi)+xs, OP_wm_temp, c*x)
xlim([0.5, Nroi+0.5])
set(gca,'XTick',[],'YTick',[0 0.5 1]);
title('FA, uFA, OP')
subplot(4,1,3),axis([0.5 Nroi+0.5 0 1])
for r = 1:Nroi
    t = text(r,0,new_names{r},'Rotation',90,'FontSize',10);
end
axis off

% individual values
sub_mat = zeros(Len, Nroi);
FA_wm_mat=zeros(Len, Nroi);
MD_wm_mat=zeros(Len, Nroi);
AD_wm_mat=zeros(Len, Nroi);
RD_wm_mat=zeros(Len, Nroi);
MK_wm_mat=zeros(Len, Nroi);
AK_wm_mat=zeros(Len, Nroi);
RK_wm_mat=zeros(Len, Nroi);
MKv_wm_mat=zeros(Len, Nroi);
AKv_wm_mat=zeros(Len, Nroi);
RKv_wm_mat=zeros(Len, Nroi);
MKu_wm_mat=zeros(Len, Nroi);
AKu_wm_mat=zeros(Len, Nroi);
RKu_wm_mat=zeros(Len, Nroi);

uFA_wm_mat=zeros(Len, Nroi);
OP_wm_mat=zeros(Len, Nroi);

for s=1:Len
    [path, dwi_name, bvaln, bvecn, mn, dt] = fun_data_dir(s);
    load([path, fs, 'cti_gs_c'])
    
    V = load_untouch_nii([path, fs, fold, fs, 'ROIs.nii']);
    WM=V.img;
    for ri = 1:length(rois)
        ri_ind = rois(ri);
        if ri_ind < 7
            WMroi = WM == ri_ind;
        else
            WMroi = (WM == ri_ind) + (WM == ri_ind+1);
        end
        FA_wm_mat(s, ri) = mean(FA(WMroi>0));
        MD_wm_mat(s, ri) = mean(MD(WMroi>0));
        AD_wm_mat(s, ri) = mean(AD(WMroi>0));
        RD_wm_mat(s, ri) = mean(RD(WMroi>0));
        
        
        MKT(MKT(:)<0)=0;
        MKT(MKT(:)>3)=3;
        AKT(AKT(:)<0)=0;
        AKT(AKT(:)>3)=3;
        RKT(RKT(:)<0)=0;
        RKT(RKT(:)>3)=3;
        
        MKTv(MKTv(:)<0)=0;
        MKTv(MKTv(:)>3)=3;
        AKTv(AKTv(:)<0)=0;
        AKTv(AKTv(:)>3)=3;
        RKTv(RKTv(:)<0)=0;
        RKTv(RKTv(:)>3)=3;
        
        MKTi(MKTi(:)<0)=0;
        MKTi(MKTi(:)>3)=3;
        AKTi(AKTi(:)<0)=0;
        AKTi(AKTi(:)>3)=3;
        RKTi(RKTi(:)<0)=0;
        RKTi(RKTi(:)>3)=3;
        
        MK_wm_mat(s, ri) = mean(MKT(WMroi>0));
        AK_wm_mat(s, ri) = mean(AKT(WMroi>0));
        RK_wm_mat(s, ri) = mean(RKT(WMroi>0));
        
        MKv_wm_mat(s, ri) = mean(MKTv(WMroi>0));
        AKv_wm_mat(s, ri) = mean(AKTv(WMroi>0));
        RKv_wm_mat(s, ri) = mean(RKTv(WMroi>0));
        
        MKu_wm_mat(s, ri) = mean(MKTi(WMroi>0));
        AKu_wm_mat(s, ri) = mean(AKTi(WMroi>0));
        RKu_wm_mat(s, ri) = mean(RKTi(WMroi>0));
        
        OP = order_parameter(KANISO, FA);
        OP_wm_mat(s, ri) = mean(OP(WMroi>0));
        
        numerator = 15.0 * KANISO;
        denominator = 10.0 * KANISO + 12.0;

        ufa2 = numerator ./ denominator;

        ufa2 = max(min(ufa2, 1.0), 0.0);
        
        uFA = sqrt(ufa2);
        
        uFA_wm_mat(s, ri) = mean(uFA(WMroi>0));

        sub_mat(s, ri) = ri;
    end
end

OP_wm_temp = mean(OP_wm_mat, 1);


figure('color', [1 1 1])
subplot(2,1,1)
bar((1:Nroi)-xs, MD_wm_temp, c*x)
hold on
bar(1:Nroi, RD_wm_temp, c*x)
bar((1:Nroi)+xs, AD_wm_temp, c*x)
plot(sub_mat(:)-xs, MD_wm_mat(:), '.black')
plot(sub_mat(:), RD_wm_mat(:), '.black')
plot(sub_mat(:)+xs, AD_wm_mat(:), '.black')
xlim([0.5, Nroi+0.5])
set(gca,'XTick',[],'YTick',[0 0.5 1 1.5 2]);
title('$$D$$', 'Interpreter','latex', 'fontsize', fff)
legend('$$\overline{D}$$', '$$D_\bot$$', '$$D_\parallel$$',...
    'Interpreter','latex', 'fontsize', fff)
subplot(4,1,3),axis([0.5 Nroi+0.5 0 1])
for r = 1:Nroi
    t = text(r,0,new_names{r},'Rotation',90,'FontSize',10);
end
axis off



figure('color', [1 1 1])
subplot(2,1,1)
bar((1:Nroi)-xs, MK_wm_temp, c*x)
hold on
bar(1:Nroi, RK_wm_temp, c*x)
bar((1:Nroi)+xs, AK_wm_temp, c*x)
plot(sub_mat(:)-xs, MK_wm_mat(:), '.black')
plot(sub_mat(:), RK_wm_mat(:), '.black')
plot(sub_mat(:)+xs, AK_wm_mat(:), '.black')
xlim([0.5, Nroi+0.5])
set(gca,'XTick',[],'YTick',[0 0.5 1 1.5 2]);
title('$$W$$', 'Interpreter','latex', 'fontsize', fff)
legend('$$\overline{W}$$', '$$W_\bot$$', '$$W_\parallel$$',...
    'Interpreter','latex', 'fontsize', fff)
subplot(4,1,3),axis([0.5 Nroi+0.5 0 1])
for r = 1:Nroi
    t = text(r,0,new_names{r},'Rotation',90,'FontSize',10);
end
axis off


figure('color', [1 1 1])
subplot(2,1,1)
bar((1:Nroi)-xs, MKv_wm_temp, c*x)
hold on
bar(1:Nroi, RKv_wm_temp, c*x)
bar((1:Nroi)+xs, AKv_wm_temp, c*x)
plot(sub_mat(:)-xs, MKv_wm_mat(:), '.black')
plot(sub_mat(:), RKv_wm_mat(:), '.black')
plot(sub_mat(:)+xs, AKv_wm_mat(:), '.black')
xlim([0.5, Nroi+0.5])
set(gca,'XTick',[],'YTick',[0 0.5 1 1.5 2]);
title('$$vW$$', 'Interpreter','latex', 'fontsize', fff)
legend('$$v\overline{W}$$', '$$vW_\bot$$', '$$vW_\parallel$$',...
    'Interpreter','latex', 'fontsize', fff)
subplot(4,1,3),axis([0.5 Nroi+0.5 0 1])
for r = 1:Nroi
    t = text(r,0,new_names{r},'Rotation',90,'FontSize',10);
end
axis off

figure('color', [1 1 1])
subplot(2,1,1)
bar((1:Nroi)-xs, MKu_wm_temp, c*x)
hold on
bar(1:Nroi, RKu_wm_temp, c*x)
bar((1:Nroi)+xs, AKu_wm_temp, c*x)
plot(sub_mat(:)-xs, MKu_wm_mat(:), '.black')
plot(sub_mat(:), RKu_wm_mat(:), '.black')
plot(sub_mat(:)+xs, AKu_wm_mat(:), '.black')
xlim([0.5, Nroi+0.5])
set(gca,'XTick',[],'YTick',[0 0.5 1 1.5 2]);
title('$${\mu}W$$', 'Interpreter','latex', 'fontsize', fff)
legend('$${\mu}\overline{W}$$', '$${\mu}W_\bot$$', '$${\mu}W_\parallel$$',...
    'Interpreter','latex', 'fontsize', fff)
subplot(4,1,3),axis([0.5 Nroi+0.5 0 1])
for r = 1:Nroi
    t = text(r,0,new_names{r},'Rotation',90,'FontSize',10);
end
axis off

figure('color', [1 1 1])
subplot(2,1,1)
bar((1:Nroi)-xs, FA_wm_temp, c*x)
hold on
bar(1:Nroi, uFA_wm_temp, c*x)
bar((1:Nroi)+xs, OP_wm_temp, c*x)
plot(sub_mat(:)-xs, FA_wm_mat(:), '.black')
plot(sub_mat(:), uFA_wm_mat(:), '.black')
plot(sub_mat(:)+xs, OP_wm_mat(:), '.black')
xlim([0.5, Nroi+0.5])
set(gca,'XTick',[],'YTick',[0 0.5 1 1.5]);
title('FA/uFA/OP', 'fontsize', fff)
legend('$${\mu}\overline{W}$$', '$${\mu}W_\bot$$', '$${\mu}W_\parallel$$',...
    'Interpreter','latex', 'fontsize', fff)
subplot(4,1,3),axis([0.5 Nroi+0.5 0 1])
for r = 1:Nroi
    t = text(r,0,new_names{r},'Rotation',90,'FontSize',10);
end
axis off

close all

%% sort rois
if sort_M == 1
    [m, ind] = sort(mean(OP_wm_mat,1));
else
    [m, ind] = sort(RKu_wm_temp);
end
rois = rois(ind);

pth = 'C:\Users\rafae\Data\CTIhuman\sub-';

nvoxel_wm_mat = zeros(Len, length(rois));
for s = Len:-1:1
    if s < 10
        V = load_untouch_nii([pth, '0', num2str(s), fs, fold, fs, 'ROIs.nii']);
    else
        V = load_untouch_nii([pth, num2str(s), fs, fold, fs, 'ROIs.nii']);
    end
    WM=V.img;
    for ri = 1:length(rois)
        ri_ind = rois(ri);
        if ri_ind < 7
            nvoxel_wm_mat(s, ri) = sum(WM(:) == ri_ind);
        else
            nvoxel_wm_mat(s, ri) = sum(WM(:) == ri_ind) + sum(WM(:) == ri_ind+1);
        end
    end
end


% Prepare names
Nroi = 0; new_names = {};
for r = rois % unimodal ROIs (i.e. 6 first ROIs)
    Nroi = Nroi+1;
    if r < 7
        name = roi_names{r};
        name(find(name=='_')) = ' ';
    else
        name = roi_names{r}(1:(end-2));
    end
    name(find(name == '_')) = ' ';
    new_names{Nroi} = name;
end

% average number of voxels
mnv = mean(nvoxel_wm_mat, 1);
figure('color', [1 1 1]);
subplot(2,1,1),axis([0.5 27.5 0 1])
plot(mnv)
yline(100)
xlim([0.5, 27.5])
%set(gca, 'Xtick', 1:27, 'XTickLabels', {new_names})
subplot(2,1,2),axis([0.5 27.5 0 1])
for r = 1:Nroi
    t = text(r,0,new_names{r},'Rotation',90,'FontSize',10);
end
axis off

% template values
FA_wm_temp=zeros(1, Nroi);
MD_wm_temp=zeros(1, Nroi);
AD_wm_temp=zeros(1, Nroi);
RD_wm_temp=zeros(1, Nroi);
MK_wm_temp=zeros(1, Nroi);
AK_wm_temp=zeros(1, Nroi);
RK_wm_temp=zeros(1, Nroi);
MKv_wm_temp=zeros(1, Nroi);
AKv_wm_temp=zeros(1, Nroi);
RKv_wm_temp=zeros(1, Nroi);
MKu_wm_temp=zeros(1, Nroi);
AKu_wm_temp=zeros(1, Nroi);
RKu_wm_temp=zeros(1, Nroi);
uFA_wm_temp=zeros(1, Nroi);

load_templated

for ri = 1:length(rois)
    ri_ind = rois(ri);
    if ri_ind < 7
        WMroi = ROI_temp(:) == ri_ind;
    else
        WMroi = (ROI_temp(:) == ri_ind) + (ROI_temp(:) == ri_ind+1);
    end
    FA_wm_temp(ri) = mean(FA(WMroi>0));
    MD_wm_temp(ri) = mean(MD(WMroi>0));
    AD_wm_temp(ri) = mean(AD(WMroi>0));
    RD_wm_temp(ri) = mean(RD(WMroi>0));
    
    MK_wm_temp(ri) = mean(MK(WMroi>0));
    AK_wm_temp(ri) = mean(AK(WMroi>0));
    RK_wm_temp(ri) = mean(RK(WMroi>0));
    
    MKv_wm_temp(ri) = mean(MKv(WMroi>0));
    AKv_wm_temp(ri) = mean(AKv(WMroi>0));
    RKv_wm_temp(ri) = mean(RKv(WMroi>0));
    
    MKu_wm_temp(ri) = mean(MKi(WMroi>0));
    AKu_wm_temp(ri) = mean(AKi(WMroi>0));
    RKu_wm_temp(ri) = mean(RKi(WMroi>0));
    
    uFA_wm_temp(ri) = mean(uFA(WMroi>0));
    
end
c = 4;
x = 1/(3*c+6);
xs = (c+1)*x;
figure('color', [1 1 1])
subplot(2,1,1)
bar((1:Nroi)-xs, MD_wm_temp, c*x)
hold on
bar(1:Nroi, RD_wm_temp, c*x)
bar((1:Nroi)+xs, AD_wm_temp, c*x)
xlim([0.5, Nroi+0.5])
set(gca,'XTick',[],'YTick',[0 0.5 1 1.5 2]);
title('Diff')
subplot(4,1,3),axis([0.5 Nroi+0.5 0 1])
for r = 1:Nroi
    t = text(r,0,new_names{r},'Rotation',90,'FontSize',10);
end
axis off


figure('color', [1 1 1])
subplot(2,1,1)
bar((1:Nroi)-xs, MK_wm_temp, c*x)
hold on
bar(1:Nroi, RK_wm_temp, c*x)
bar((1:Nroi)+xs, AK_wm_temp, c*x)
xlim([0.5, Nroi+0.5])
set(gca,'XTick',[],'YTick',[0 0.5 1 1.5 2]);
title('Kurt')
subplot(4,1,3),axis([0.5 Nroi+0.5 0 1])
for r = 1:Nroi
    t = text(r,0,new_names{r},'Rotation',90,'FontSize',10);
end
axis off


figure('color', [1 1 1])
subplot(2,1,1)
bar((1:Nroi)-xs, MKv_wm_temp, c*x)
hold on
bar(1:Nroi, RKv_wm_temp, c*x)
bar((1:Nroi)+xs, AKv_wm_temp, c*x)
xlim([0.5, Nroi+0.5])
set(gca,'XTick',[],'YTick',[0 0.5 1 1.5 2]);
title('Var Kurt')
subplot(4,1,3),axis([0.5 Nroi+0.5 0 1])
for r = 1:Nroi
    t = text(r,0,new_names{r},'Rotation',90,'FontSize',10);
end
axis off

figure('color', [1 1 1])
subplot(2,1,1)
bar((1:Nroi)-xs, MKu_wm_temp, c*x)
hold on
bar(1:Nroi, RKu_wm_temp, c*x)
bar((1:Nroi)+xs, AKu_wm_temp, c*x)
xlim([0.5, Nroi+0.5])
set(gca,'XTick',[],'YTick',[0 0.5 1 1.5 2]);
title('Micro Kurt')
subplot(4,1,3),axis([0.5 Nroi+0.5 0 1])
for r = 1:Nroi
    t = text(r,0,new_names{r},'Rotation',90,'FontSize',10);
end
axis off


figure('color', [1 1 1])
subplot(2,1,1)
bar((1:Nroi)-xs, MKv_wm_temp./MK_wm_temp, c*x)
hold on
bar(1:Nroi, RKv_wm_temp./RK_wm_temp, c*x)
bar((1:Nroi)+xs, AKv_wm_temp./AK_wm_temp, c*x)
xlim([0.5, Nroi+0.5])
set(gca,'XTick',[],'YTick',[0 0.5 1]);
title('Kv/Kt')
subplot(4,1,3),axis([0.5 Nroi+0.5 0 1])
for r = 1:Nroi
    t = text(r,0,new_names{r},'Rotation',90,'FontSize',10);
end
axis off

% individual values
sub_mat = zeros(Len, Nroi);
FA_wm_mat=zeros(Len, Nroi);
MD_wm_mat=zeros(Len, Nroi);
AD_wm_mat=zeros(Len, Nroi);
RD_wm_mat=zeros(Len, Nroi);
MK_wm_mat=zeros(Len, Nroi);
AK_wm_mat=zeros(Len, Nroi);
RK_wm_mat=zeros(Len, Nroi);
MKv_wm_mat=zeros(Len, Nroi);
AKv_wm_mat=zeros(Len, Nroi);
RKv_wm_mat=zeros(Len, Nroi);
MKu_wm_mat=zeros(Len, Nroi);
AKu_wm_mat=zeros(Len, Nroi);
RKu_wm_mat=zeros(Len, Nroi);

uFA_wm_mat=zeros(Len, Nroi);
OP_wm_mat=zeros(Len, Nroi);

for s=1:Len
    [path, dwi_name, bvaln, bvecn, mn, dt] = fun_data_dir(s);
    load([path, fs, 'cti_gs_c'])
    
    
    
    V = load_untouch_nii([path, fs, fold, fs, 'ROIs.nii']);
    WM=V.img;
    for ri = 1:length(rois)
        ri_ind = rois(ri);
        if ri_ind < 7
            WMroi = WM == ri_ind;
        else
            WMroi = (WM == ri_ind) + (WM == ri_ind+1);
        end
        FA_wm_mat(s, ri) = mean(FA(WMroi>0));
        MD_wm_mat(s, ri) = mean(MD(WMroi>0));
        AD_wm_mat(s, ri) = mean(AD(WMroi>0));
        RD_wm_mat(s, ri) = mean(RD(WMroi>0));
        
        
        MKT(MKT(:)<0)=0;
        MKT(MKT(:)>3)=3;
        AKT(AKT(:)<0)=0;
        AKT(AKT(:)>3)=3;
        RKT(RKT(:)<0)=0;
        RKT(RKT(:)>3)=3;
        
        MKTv(MKTv(:)<0)=0;
        MKTv(MKTv(:)>3)=3;
        AKTv(AKTv(:)<0)=0;
        AKTv(AKTv(:)>3)=3;
        RKTv(RKTv(:)<0)=0;
        RKTv(RKTv(:)>3)=3;
        
        MKTi(MKTi(:)<0)=0;
        MKTi(MKTi(:)>3)=3;
        AKTi(AKTi(:)<0)=0;
        AKTi(AKTi(:)>3)=3;
        RKTi(RKTi(:)<0)=0;
        RKTi(RKTi(:)>3)=3;
        
        MK_wm_mat(s, ri) = mean(MKT(WMroi>0));
        AK_wm_mat(s, ri) = mean(AKT(WMroi>0));
        RK_wm_mat(s, ri) = mean(RKT(WMroi>0));
        
        MKv_wm_mat(s, ri) = mean(MKTv(WMroi>0));
        AKv_wm_mat(s, ri) = mean(AKTv(WMroi>0));
        RKv_wm_mat(s, ri) = mean(RKTv(WMroi>0));
        
        MKu_wm_mat(s, ri) = mean(MKTi(WMroi>0));
        AKu_wm_mat(s, ri) = mean(AKTi(WMroi>0));
        RKu_wm_mat(s, ri) = mean(RKTi(WMroi>0));
        
        OP = order_parameter(KANISO, FA);
        OP_wm_mat(s, ri) = mean(OP(WMroi>0));
        
        numerator = 15.0 * KANISO;
        denominator = 10.0 * KANISO + 12.0;

        ufa2 = numerator ./ denominator;

        ufa2 = max(min(ufa2, 1.0), 0.0);
        
        uFA = sqrt(ufa2);
        
        uFA_wm_mat(s, ri) = mean(uFA(WMroi>0));
        
        sub_mat(s, ri) = ri;
    end
end

OP_wm_temp = mean(OP_wm_mat, 1);

c = 6;
x = 1/(2*c+5);
xs = (c+1)*x/2;

figure('color', [1 1 1])
subplot(2,1,1)
bar((1:Nroi)-xs, RD_wm_temp, c*x)
hold on
bar((1:Nroi)+xs, AD_wm_temp, c*x)
plot(sub_mat(:)-xs, RD_wm_mat(:), '.black')
plot(sub_mat(:)+xs, AD_wm_mat(:), '.black')
xlim([0.5, Nroi+0.5])
set(gca,'XTick',[],'YTick',[0 0.5 1 1.5 2]);
title('$$D$$', 'Interpreter','latex', 'fontsize', fff)
legend('$$D_\bot$$', '$$D_\parallel$$',...
    'Interpreter','latex', 'fontsize', fff)
subplot(4,1,3),axis([0.5 Nroi+0.5 0 1])
for r = 1:Nroi
    t = text(r,0,new_names{r},'Rotation',90,'FontSize',10);
end
axis off



figure('color', [1 1 1])
subplot(2,1,1)
bar((1:Nroi)-xs, RK_wm_temp, c*x)
hold on
bar((1:Nroi)+xs, AK_wm_temp, c*x)
plot(sub_mat(:)-xs, RK_wm_mat(:), '.black')
plot(sub_mat(:)+xs, AK_wm_mat(:), '.black')
xlim([0.5, Nroi+0.5])
set(gca,'XTick',[],'YTick',[0 0.5 1 1.5 2]);
%title('$$W$$', 'Interpreter','latex', 'fontsize', fff)
legend('$$W_\bot$$', '$$W_\parallel$$',...
    'Interpreter','latex', 'fontsize', fff)
subplot(4,1,3),axis([0.5 Nroi+0.5 0 1])
for r = 1:Nroi
    t = text(r,0,new_names{r},'Rotation',90,'FontSize',10);
end
axis off

figure('color', [1 1 1])
subplot(2,1,1)
bar((1:Nroi)-xs, FA_wm_temp, c*x)
hold on
bar(1:Nroi, uFA_wm_temp, c*x)
bar((1:Nroi)+xs, OP_wm_temp, c*x)
plot(sub_mat(:)-xs, FA_wm_mat(:), '.black')
plot(sub_mat(:), uFA_wm_mat(:), '.black')
plot(sub_mat(:)+xs, OP_wm_mat(:), '.black')
xlim([0.5, Nroi+0.5])
set(gca,'XTick',[],'YTick',[0 0.5 1 1.5]);
title('FA/uFA/OP', 'fontsize', fff)
legend('$${\mu}\overline{W}$$', '$${\mu}W_\bot$$', '$${\mu}W_\parallel$$',...
    'Interpreter','latex', 'fontsize', fff)
subplot(4,1,3),axis([0.5 Nroi+0.5 0 1])
for r = 1:Nroi
    t = text(r,0,new_names{r},'Rotation',90,'FontSize',10);
end
axis off


figure('color', [1 1 1])
subplot(2,1,1)
bar((1:Nroi)-xs, RKv_wm_temp, c*x)
hold on
bar((1:Nroi)+xs, AKv_wm_temp, c*x)
plot(sub_mat(:)-xs, RKv_wm_mat(:), '.black')
plot(sub_mat(:)+xs, AKv_wm_mat(:), '.black')
xlim([0.5, Nroi+0.5])
ylim([0 2.5])
set(gca,'XTick',[],'YTick',[0 0.5 1 1.5 2]);
%title('$$vW$$', 'Interpreter','latex', 'fontsize', fff)
legend('$$vW_\bot$$', '$$vW_\parallel$$',...
    'Interpreter','latex', 'fontsize', fff)
subplot(4,1,3),axis([0.5 Nroi+0.5 0 1])
for r = 1:Nroi
    t = text(r,0,new_names{r},'Rotation',90,'FontSize',10);
end
axis off

for ri = 1:19
    deltaKv = RKv_wm_mat(:, ri) - AKv_wm_mat(:, ri);
    [h,pval,ci,stats] = ttest(deltaKv);
    pvar(ri) = pval;
end
[Q, ~, ~, adj_p] = fdr_bh(pvar(:));

figure('color', [1 1 1])
subplot(2,1,1)
bar((1:Nroi)-xs, RKu_wm_temp, c*x)
hold on
bar((1:Nroi)+xs, AKu_wm_temp, c*x)
plot(sub_mat(:)-xs, RKu_wm_mat(:), '.black')
plot(sub_mat(:)+xs, AKu_wm_mat(:), '.black')
xlim([0.5, Nroi+0.5])
set(gca,'XTick',[],'YTick',[0 0.1 0.2 0.3 0.4 0.5]);
%title('$${\mu}W$$', 'Interpreter','latex', 'fontsize', fff)
legend('$${\mu}W_\bot$$', '$${\mu}W_\parallel$$',...
    'Interpreter','latex', 'fontsize', fff)
subplot(4,1,3),axis([0.5 Nroi+0.5 0 1])
for r = 1:Nroi
    t = text(r,0,new_names{r},'Rotation',90,'FontSize',10);
end
axis off

for ri = 1:19
    deltaKu = RKu_wm_mat(:, ri) - AKu_wm_mat(:, ri);
    [h,pval,ci,stats] = ttest(deltaKu);
    pmic(ri) = pval;
end
[Q, ~, ~, adj_p] = fdr_bh(pmic(:));

figure('color', [1 1 1])
subplot(2,1,1)
bar((1:Nroi)-xs, RKu_wm_temp./RK_wm_temp, c*x)
hold on
bar((1:Nroi)+xs, AKu_wm_temp./AK_wm_temp, c*x)
plot(sub_mat(:)-xs, RKu_wm_mat(:)./RK_wm_mat(:), '.black')
plot(sub_mat(:)+xs, AKu_wm_mat(:)./AK_wm_mat(:), '.black')
xlim([0.5, Nroi+0.5])
set(gca,'XTick',[],'YTick',[0 0.05 0.1 0.15 0.2 0.25]);
%title('$${\mu}W/W$$', 'Interpreter','latex', 'fontsize', fff)
legend('$${\mu}W_\bot/W_\bot$$', '$${\mu}W_\parallel/W_\parallel$$',...
    'Interpreter','latex', 'fontsize', fff, 'Location', 'northwest')
subplot(4,1,3),axis([0.5 Nroi+0.5 0 1])
for r = 1:Nroi
    t = text(r,0,new_names{r},'Rotation',90,'FontSize',10);
end
axis off


for ri = 1:19
    deltaKu = RKu_wm_mat(:, ri)./RK_wm_mat(:,ri) - AKu_wm_mat(:, ri)./AK_wm_mat(:,ri);
    [h,pval,ci,stats] = ttest(deltaKu);
    pmic(ri) = pval;
end
[Q, ~, ~, adj_p] = fdr_bh(pmic(:));

figure('color', [1 1 1])
subplot(2,1,1)
bar((1:Nroi)-xs, RKv_wm_temp./RK_wm_temp, c*x)
hold on
bar((1:Nroi)+xs, AKv_wm_temp./AK_wm_temp, c*x)
plot(sub_mat(:)-xs, RKv_wm_mat(:)./RK_wm_mat(:), '.black')
plot(sub_mat(:)+xs, AKv_wm_mat(:)./AK_wm_mat(:), '.black')
xlim([0.5, Nroi+0.5])
set(gca,'XTick',[],'YTick',[0 0.25 0.5 0.75 1]);
%title('$$vW/W$$', 'Interpreter','latex', 'fontsize', fff)
legend('$$vW_\bot/W_\bot$$', '$$vW_\parallel/W_\parallel$$',...
    'Interpreter','latex', 'fontsize', fff, 'Location', 'northwest')
subplot(4,1,3),axis([0.5 Nroi+0.5 0 1])
for r = 1:Nroi
    t = text(r,0,new_names{r},'Rotation',90,'FontSize',10);
end
axis off

save('roi_analysis2', 'new_names', 'Nroi', 'xs', 'sub_mat', 'c', ...
    'FA_wm_temp', 'FA_wm_mat',...
    'RKv_wm_temp', 'RKu_wm_temp', 'RK_wm_temp', ...
    'AKv_wm_temp', 'AKu_wm_temp', 'AK_wm_temp', 'OP_wm_temp',...
    'RKv_wm_mat', 'RKu_wm_mat', 'RK_wm_mat', ...
    'AKv_wm_mat', 'AKu_wm_mat', 'AK_wm_mat', 'OP_wm_mat')