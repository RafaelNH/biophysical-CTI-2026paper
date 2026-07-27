addpath(['..',fs,'DWIu_toolbox_v1.11',fs,'DDE_toolbox'])

WD = 1:11;
ns = length(WD);

% remove influence of negative kurtosis values than can substancially
% leverage the analysis?
rnegk = true;

%% Varable Initialization
% diffusion quantities initialization
fa = zeros(1, ns);
md = zeros(1, ns);
rd = zeros(1, ns);
ad = zeros(1, ns);

% total kurtosis quantities initialization
mk = zeros(1, ns);
rk = zeros(1, ns);
ak = zeros(1, ns);

% variance kurtosis quantities initialization
vk = zeros(1, ns);
vrk = zeros(1, ns);
vak = zeros(1, ns);

% microscopic kurtosis quantities initialization
uk = zeros(1, ns);
urk = zeros(1, ns);
uak = zeros(1, ns);

% CTI Kaniso/Kiso/Css
kaniso = zeros(1, ns);
kiso = zeros(1, ns);
cov = zeros(1, ns);

% micro-anisotropic/dispersion quantities
op = zeros(1, ns);
ufa = zeros(1, ns);

%sss = [6, 6, 7, 7, 6, 6, 7, 6, 6, 6, 6];
controls = [2, 5, 6, 8, 9, 11];
tumor = [1, 3, 4, 7, 10];
control_logic = false(1, 11);
control_logic(controls) = true;

for ti = WD
    
    dinfo = fun_data_dirs(ti);
    load([dinfo(1).savename(1:end-6), 'cti__complex_denalign_drift_g_c'])
    load([dinfo(1).savename(1:end-6), 'MASK_cerebellum.mat'], 'mask_cerebellum')
    
    % Calculating OP and uFA
    numerator = 15.0 * KANISO;
    denominator = 10.0 * KANISO + 12.0;

    ufa2 = numerator ./ denominator;
    ufa2 = max(min(ufa2, 1.0), 0.0);
    uFA = sqrt(ufa2);
    
    OP = order_parameter_from_ufa(uFA, FA);
    OP(FA>uFA) = 0;
    
    
    %% remove the influence of very negative kurtosis values (if selected)
    if rnegk
        mask_cerebellum(MKT(:)<0)=0;
        mask_cerebellum(MKT(:)>3)=0;
        mask_cerebellum(AKT(:)<0)=0;
        mask_cerebellum(AKT(:)>3)=0;
        mask_cerebellum(RKT(:)<0)=0;
        mask_cerebellum(RKT(:)>3)=0;
        
        mask_cerebellum(MKTv(:)<0)=0;
        mask_cerebellum(MKTv(:)>3)=0;
        mask_cerebellum(AKTv(:)<0)=0;
        mask_cerebellum(AKTv(:)>3)=0;
        mask_cerebellum(RKTv(:)<0)=0;
        mask_cerebellum(RKTv(:)>3)=0;
        
        mask_cerebellum(MKTi(:)<0)=0;
        mask_cerebellum(MKTi(:)>3)=0;
        mask_cerebellum(AKTi(:)<0)=0;
        mask_cerebellum(AKTi(:)>3)=0;
        mask_cerebellum(RKTi(:)<0)=0;
        mask_cerebellum(RKTi(:)>3)=0;
        
    end
    
    relAKTi = AKTi./AKT;
    relRKTi = RKTi./RKT;
    relAKTv = AKTv./AKT;
    relRKTv = RKTv./RKT;
    
    if ti == 2
        si = 5;
        FAcontrol = FA(end:-1:1, :, si);
        MDcontrol = MD(end:-1:1, :, si);
        RDcontrol = RD(end:-1:1, :, si);
        ADcontrol = AD(end:-1:1, :, si);
        
        MKcontrol = MKT(end:-1:1, :, si);
        RKcontrol = RKT(end:-1:1, :, si);
        AKcontrol = AKT(end:-1:1, :, si);
        
        MKicontrol = MKTi(end:-1:1, :, si);
        RKicontrol = RKTi(end:-1:1, :, si);
        AKicontrol = AKTi(end:-1:1, :, si);
        
        MKvcontrol = MKTv(end:-1:1, :, si);
        RKvcontrol = RKTv(end:-1:1, :, si);
        AKvcontrol = AKTv(end:-1:1, :, si);
        
        KANIcontrol = KANISO(end:-1:1, :, si);
        KISOcontrol = KISO(end:-1:1, :, si);
        COVcontrol = Cov_ss_n(end:-1:1, :, si);
        
        uFAcontrol = uFA(end:-1:1, :, si);
        OPcontrol = OP(end:-1:1, :, si);
        
        relAKTicontrol = relAKTi(end:-1:1, :, si);
        relRKTicontrol = relRKTi(end:-1:1, :, si);
        
        relAKTvcontrol = relAKTv(end:-1:1, :, si);
        relRKTvcontrol = relRKTv(end:-1:1, :, si);
        
    elseif ti==1
        si = 6;
        FAtumor = FA(end:-1:1, :, si);
        MDtumor = MD(end:-1:1, :, si);
        RDtumor = RD(end:-1:1, :, si);
        ADtumor = AD(end:-1:1, :, si);
        
        MKtumor = MKT(end:-1:1, :, si);
        RKtumor = RKT(end:-1:1, :, si);
        AKtumor = AKT(end:-1:1, :, si);
        
        MKitumor = MKTi(end:-1:1, :, si);
        RKitumor = RKTi(end:-1:1, :, si);
        AKitumor = AKTi(end:-1:1, :, si);
        
        MKvtumor = MKTv(end:-1:1, :, si);
        RKvtumor = RKTv(end:-1:1, :, si);
        AKvtumor = AKTv(end:-1:1, :, si);
        
        KANItumor = KANISO(end:-1:1, :, si);
        KISOtumor = KISO(end:-1:1, :, si);
        COVtumor = Cov_ss_n(end:-1:1, :, si);
        
        relAKTitumor = relAKTi(end:-1:1, :, si);
        relRKTitumor = relRKTi(end:-1:1, :, si);
        
        relAKTvtumor = relAKTv(end:-1:1, :, si);
        relRKTvtumor = relRKTv(end:-1:1, :, si);
        
        uFAtumor = uFA(end:-1:1, :, si);
        OPtumor = OP(end:-1:1, :, si);
    end
    
    FAv = FA(mask_cerebellum(:)==1);
    MDv = MD(mask_cerebellum(:)==1);
    RDv = RD(mask_cerebellum(:)==1);
    ADv = AD(mask_cerebellum(:)==1);
    
    MKv = MKT(mask_cerebellum(:)==1);
    RKv = RKT(mask_cerebellum(:)==1);
    AKv = AKT(mask_cerebellum(:)==1);
    
    MKTvv = MKTv(mask_cerebellum(:)==1);
    RKTvv = RKTv(mask_cerebellum(:)==1);
    AKTvv = AKTv(mask_cerebellum(:)==1);
    
    MKTiv = MKTi(mask_cerebellum(:)==1);
    RKTiv = RKTi(mask_cerebellum(:)==1);
    AKTiv = AKTi(mask_cerebellum(:)==1);
    
    KANISOv = KANISO(mask_cerebellum(:)==1);
    KISOv = KISO(mask_cerebellum(:)==1);
    COVv = Cov_ss_n(mask_cerebellum(:)==1);
    
    uFAv = uFA(mask_cerebellum(:)==1);
    OPv = OP(mask_cerebellum(:)==1);
    
    relRKTvv = relRKTv(mask_cerebellum(:)==1);
    relAKTvv = relAKTv(mask_cerebellum(:)==1);
    relRKTiv = relRKTi(mask_cerebellum(:)==1);
    relAKTiv = relAKTi(mask_cerebellum(:)==1);
    
    disp(sum(mask_cerebellum(:)))
    
    % for histograms
    cerebellums(ti).FAv = FAv;
    cerebellums(ti).MDv = MDv;
    cerebellums(ti).RDv = RDv;
    cerebellums(ti).ADv = ADv;
    
    cerebellums(ti).MKv = MKv;
    cerebellums(ti).RKv = RKv;
    cerebellums(ti).AKv = AKv;
    
    cerebellums(ti).MKTiv = MKTiv;
    cerebellums(ti).RKTiv = RKTiv;
    cerebellums(ti).AKTiv = AKTiv;
    
    cerebellums(ti).MKTvv = MKTvv;
    cerebellums(ti).RKTvv = RKTvv;
    cerebellums(ti).AKTvv = AKTvv;
    
    cerebellums(ti).KANISOv = KANISOv;
    cerebellums(ti).KISOv = KISOv;
    cerebellums(ti).COVv = COVv;
    
    cerebellums(ti).uFAv = uFAv;
    cerebellums(ti).OPv = OPv;
    
    cerebellums(ti).relRKTvv = relRKTvv;
    cerebellums(ti).relAKTvv = relAKTvv;
    cerebellums(ti).relRKTiv = relRKTiv;
    cerebellums(ti).relAKTiv = relAKTiv;
    
    % for bar plots
    % Diffusion quantities
    fa(ti) = mean(FA(mask_cerebellum(:)==1));
    md(ti) = mean(MD(mask_cerebellum(:)==1));
    rd(ti) = mean(RD(mask_cerebellum(:)==1));
    ad(ti) = mean(AD(mask_cerebellum(:)==1));
    
    % total kurtosis quantities
    mk(ti) = mean(MKT(mask_cerebellum(:)==1));
    rk(ti) = mean(RKT(mask_cerebellum(:)==1));
    ak(ti) = mean(AKT(mask_cerebellum(:)==1));
    
    % variance kurtosis quantities
    vk(ti) = mean(MKTv(mask_cerebellum(:)==1));
    vrk(ti) = mean(RKTv(mask_cerebellum(:)==1));
    vak(ti) = mean(AKTv(mask_cerebellum(:)==1));
    
    % microscopic kurtosis quantities
    uk(ti) = mean(KINTRA(mask_cerebellum(:)==1));
    urk(ti) = mean(RKTi(mask_cerebellum(:)==1));
    uak(ti) = mean(AKTi(mask_cerebellum(:)==1));
    
    % CTI variance quantities kaniso/Kiso/Css
    kaniso(ti) = mean(KANISO(mask_cerebellum(:)==1));
    kiso(ti) = mean(KISO(mask_cerebellum(:)==1));
    cov(ti) = mean(Cov_ss_n(mask_cerebellum(:)==1));
    
    % Calculating OP and uFA
    op(ti) = mean(OP(mask_cerebellum(:)==1));
    ufa(ti) = mean(uFA(mask_cerebellum(:)==1));
    
end