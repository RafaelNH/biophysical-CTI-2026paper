addpath(['..',fs,'DWIu_toolbox_v1.11',fs,'DDE_toolbox'])

WD = 1:6;
ns = length(WD);
tumors = [1, 4, 6, 5, 3, 2];

% remove influence of negative kurtosis values than can substancially
% leverage the analysis?
rnegk = true;

constraints = true;
Smooth_Gaussian = false;
exlude_b0s = true;
denoised_data = true;

if constraints
    cst = '_c';
else
    cst = '';
end

if Smooth_Gaussian
    gst = '_gs';
else
    gst = '';
end

if denoised_data == true
    den = '_den';
else
    den = '';
end


ssst = 1.5;
sss = 1.5;

for ti = WD
    
    % load data
    whichd = tumors(ti);
    dinfo = fun_data_dirs(whichd);
    
    if whichd == 1 % CT1
        si = 7;
    elseif whichd == 7 % CT1
        si = 6;
    elseif whichd == 2 % GL1
        si = 14;
    elseif whichd == 3 % GL2
        si = 7;
    elseif whichd == 4 % CT2
        si = 8;
    elseif whichd == 5 % GL3
        si = 10;
    elseif whichd == 6 % CT3
        si = 10;
    end
    
    
    di = 1;
    file_name = dinfo(di).savename;
    disp(file_name)
    if whichd == 1
        load([dinfo(1).savename(1:end-7), 'cti', den, cst])
        load([dinfo(1).savename(1:end-7), 'tumour_mask'])
    elseif whichd == 4
        load([dinfo(1).savename(1:end-5), 'cti', den, cst])
        load([dinfo(1).savename(1:end-5), 'tumour_mask'])
    else
        load([dinfo(1).savename(1:end-6), 'cti', den, cst])
        load([dinfo(1).savename(1:end-6), 'tumour_mask'])
    end
    
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
        tumour_mask(MKT(:)<-0.5)=0;
        tumour_mask(MKT(:)>3)=0;
        tumour_mask(AKT(:)<-0.5)=0;
        tumour_mask(AKT(:)>3)=0;
        tumour_mask(RKT(:)<-0.5)=0;
        tumour_mask(RKT(:)>3)=0;
        
        tumour_mask(MKTv(:)<-0.5)=0;
        tumour_mask(MKTv(:)>3)=0;
        tumour_mask(AKTv(:)<-0.5)=0;
        tumour_mask(AKTv(:)>3)=0;
        tumour_mask(RKTv(:)<-0.5)=0;
        tumour_mask(RKTv(:)>3)=0;
        
        tumour_mask(MKTi(:)<-0.5)=0;
        tumour_mask(MKTi(:)>3)=0;
        tumour_mask(AKTi(:)<-0.5)=0;
        tumour_mask(AKTi(:)>3)=0;
        tumour_mask(RKTi(:)<-0.5)=0;
        tumour_mask(RKTi(:)>3)=0;
        
    end
    
    if ti == 2
        FAct2a = FA(end:-1:1, :, si);
        MDct2a = MD(end:-1:1, :, si);
        RDct2a = RD(end:-1:1, :, si);
        ADct2a = AD(end:-1:1, :, si);
        
        MKct2a = MKT(end:-1:1, :, si);
        RKct2a = RKT(end:-1:1, :, si);
        AKct2a = AKT(end:-1:1, :, si);
        
        MKict2a = MKTi(end:-1:1, :, si);
        RKict2a = RKTi(end:-1:1, :, si);
        AKict2a = AKTi(end:-1:1, :, si);
        
        MKvct2a = MKTv(end:-1:1, :, si);
        RKvct2a = RKTv(end:-1:1, :, si);
        AKvct2a = AKTv(end:-1:1, :, si);
        
        KANIct2a = KANISO(end:-1:1, :, si);
        KISOct2a = KISO(end:-1:1, :, si);
        COVct2a = Cov_ss_n(end:-1:1, :, si);
        
        uFAct2a = uFA(end:-1:1, :, si);
        OPct2a = OP(end:-1:1, :, si);
        
        
    elseif ti==4
        FAgl261 = FA(end:-1:1, :, si);
        MDgl261 = MD(end:-1:1, :, si);
        RDgl261 = RD(end:-1:1, :, si);
        ADgl261 = AD(end:-1:1, :, si);
        
        MKgl261 = MKT(end:-1:1, :, si);
        RKgl261 = RKT(end:-1:1, :, si);
        AKgl261 = AKT(end:-1:1, :, si);
        
        MKigl261 = MKTi(end:-1:1, :, si);
        RKigl261 = RKTi(end:-1:1, :, si);
        AKigl261 = AKTi(end:-1:1, :, si);
        
        MKvgl261 = MKTv(end:-1:1, :, si);
        RKvgl261 = RKTv(end:-1:1, :, si);
        AKvgl261 = AKTv(end:-1:1, :, si);
        
        KANIgl261 = KANISO(end:-1:1, :, si);
        KISOgl261 = KISO(end:-1:1, :, si);
        COVgl261 = Cov_ss_n(end:-1:1, :, si);
        
        uFAgl261 = uFA(end:-1:1, :, si);
        OPgl261 = OP(end:-1:1, :, si);
        
    end
    
    disp(sum(tumour_mask(:)))
    
    FAv = FA(tumour_mask(:)==1);
    MDv = MD(tumour_mask(:)==1);
    RDv = RD(tumour_mask(:)==1);
    ADv = AD(tumour_mask(:)==1);
    
    MKv = MKT(tumour_mask(:)==1);
    RKv = RKT(tumour_mask(:)==1);
    AKv = AKT(tumour_mask(:)==1);
    
    MKTvv = MKTv(tumour_mask(:)==1);
    RKTvv = RKTv(tumour_mask(:)==1);
    AKTvv = AKTv(tumour_mask(:)==1);
    
    MKTiv = MKTi(tumour_mask(:)==1);
    RKTiv = RKTi(tumour_mask(:)==1);
    AKTiv = AKTi(tumour_mask(:)==1);
    
    KANISOv = KANISO(tumour_mask(:)==1);
    KISOv = KISO(tumour_mask(:)==1);
    COVv = Cov_ss_n(tumour_mask(:)==1);
    
    uFAv = uFA(tumour_mask(:)==1);
    OPv = OP(tumour_mask(:)==1);
    
    gliomas(ti).MDv = MDv;
    gliomas(ti).FAv = FAv;
    gliomas(ti).RDv = RDv;
    gliomas(ti).ADv = ADv;
    
    gliomas(ti).MKv = MKv;
    gliomas(ti).RKv = RKv;
    gliomas(ti).AKv = AKv;
    
    gliomas(ti).MKTvv = MKTvv;
    gliomas(ti).RKTvv = RKTvv;
    gliomas(ti).AKTvv = AKTvv;
    
    gliomas(ti).MKTiv = MKTiv;
    gliomas(ti).RKTiv = RKTiv;
    gliomas(ti).AKTiv = AKTiv;
    
    gliomas(ti).KANISOv = KANISOv;
    gliomas(ti).KISOv = KISOv;
    gliomas(ti).COVv  = COVv;
    
    gliomas(ti).uFAv = uFAv;
    gliomas(ti).OPv = OPv;
    
    
end


