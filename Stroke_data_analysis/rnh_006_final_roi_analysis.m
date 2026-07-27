fs = filesep;
addpath(['..',fs,'DWIu_toolbox_v1.11',fs,'DDE_toolbox'])

WD = 1:10;
ns = length(WD);

% remove influence of negative kurtosis values than can substancially
% leverage the analysis?
rnegk = true;

%% Varable Initialization
% diffusion quantities initialization
fa_wm = zeros(1, ns);
fa_st_wm = zeros(1, ns);
md_wm = zeros(1, ns);
md_st_wm = zeros(1, ns);
rd_wm = zeros(1, ns);
rd_st_wm = zeros(1, ns);
ad_wm = zeros(1, ns);
ad_st_wm = zeros(1, ns);

% total kurtosis quantities initialization
mk_wm = zeros(1, ns);
mk_st_wm = zeros(1, ns);
rk_wm = zeros(1, ns);
rk_st_wm = zeros(1, ns);
ak_wm = zeros(1, ns);
ak_st_wm = zeros(1, ns);

% variance kurtosis quantities initialization
vk_wm = zeros(1, ns);
vk_st_wm = zeros(1, ns);
vrk_wm = zeros(1, ns);
vrk_st_wm = zeros(1, ns);
vak_wm = zeros(1, ns);
vak_st_wm = zeros(1, ns);

% microscopic kurtosis quantities initialization
uk_wm = zeros(1, ns);
uk_st_wm = zeros(1, ns);
urk_wm = zeros(1, ns);
urk_st_wm = zeros(1, ns);
uak_wm = zeros(1, ns);
uak_st_wm = zeros(1, ns);

% CTI Kaniso/Kiso/Css
kaniso_wm = zeros(1, ns);
kaniso_st_wm = zeros(1, ns);
kiso_wm = zeros(1, ns);
kiso_st_wm = zeros(1, ns);
cov_wm = zeros(1, ns);
cov_st_wm = zeros(1, ns);

% micro-anisotropic/dispersion quantities
op_wm = zeros(1, ns);
ufa_wm = zeros(1, ns);
op_st_wm = zeros(1, ns);
ufa_st_wm = zeros(1, ns);

%% Loading data
for ii = WD
    
    dinfo = fun_data_dirs(ii);
    if ii > 9
        load([dinfo(1).savename(1:end-6), 'dcti_d_c'])
        load([dinfo(1).savename(1:end-6), 'mask_wm_st.mat'], 'mask_wm_st')
        load([dinfo(1).savename(1:end-6), 'mask_wm_ct.mat'], 'mask_wm_ct')
        mask_wm_st(:, :, 9) = 0;
    else
        load([dinfo(1).savename(1:end-6), 'cti_d_c'])
        load([dinfo(1).savename(1:end-6), 'mask_wm_st.mat'], 'mask_wm_st')
        load([dinfo(1).savename(1:end-6), 'mask_wm_ct.mat'], 'mask_wm_ct')
    end
    
    n_wm_ct(ii) = sum(mask_wm_ct(:)==1);
    n_wm_st(ii) = sum(mask_wm_st(:)==1);
    if rnegk
        mask_wm_ct(MKT(:)<0)=0;
        mask_wm_ct(MKT(:)>3)=0;
        mask_wm_ct(AKT(:)<0)=0;
        mask_wm_ct(AKT(:)>3)=0;
        mask_wm_ct(RKT(:)<0)=0;
        mask_wm_ct(RKT(:)>3)=0;
        mask_wm_st(MKT(:)<0)=0;
        mask_wm_st(MKT(:)>3)=0;
        mask_wm_st(AKT(:)<0)=0;
        mask_wm_st(AKT(:)>3)=0;
        mask_wm_st(RKT(:)<0)=0;
        mask_wm_st(RKT(:)>3)=0;
        
        mask_wm_ct(MKTv(:)<0)=0;
        mask_wm_ct(MKTv(:)>3)=0;
        mask_wm_ct(AKTv(:)<0)=0;
        mask_wm_ct(AKTv(:)>3)=0;
        mask_wm_ct(RKTv(:)<0)=0;
        mask_wm_ct(RKTv(:)>3)=0;
        mask_wm_st(MKTv(:)<0)=0;
        mask_wm_st(MKTv(:)>3)=0;
        mask_wm_st(AKTv(:)<0)=0;
        mask_wm_st(AKTv(:)>3)=0;
        mask_wm_st(RKTv(:)<0)=0;
        mask_wm_st(RKTv(:)>3)=0;
        
        mask_wm_ct(MKTi(:)<0)=0;
        mask_wm_ct(MKTi(:)>3)=0;
        mask_wm_ct(AKTi(:)<0)=0;
        mask_wm_ct(AKTi(:)>3)=0;
        mask_wm_ct(RKTi(:)<0)=0;
        mask_wm_ct(RKTi(:)>3)=0;
        mask_wm_st(MKTi(:)<0)=0;
        mask_wm_st(MKTi(:)>3)=0;
        mask_wm_st(AKTi(:)<0)=0;
        mask_wm_st(AKTi(:)>3)=0;
        mask_wm_st(RKTi(:)<0)=0;
        mask_wm_st(RKTi(:)>3)=0;
        
    end
    n_wm_ct_f(ii) = sum(mask_wm_ct(:)==1);
    n_wm_st_f(ii) = sum(mask_wm_st(:)==1);
        
    % Calculating OP and uFA
    numerator = 15.0 * KANISO;
    denominator = 10.0 * KANISO + 12.0;

    ufa2 = numerator ./ denominator;
    ufa2 = max(min(ufa2, 1.0), 0.0);
    uFA = sqrt(ufa2);
    
    OP = order_parameter_from_ufa(uFA, FA);
    OP(FA>uFA) = 0;
    
    % Diffusion quantities
    fa_wm(ii) = mean(FA(mask_wm_ct(:)==1));
    fa_st_wm(ii) = mean(FA(mask_wm_st(:)==1));
    md_wm(ii) = mean(MD(mask_wm_ct(:)==1));
    md_st_wm(ii) = mean(MD(mask_wm_st(:)==1));
    rd_wm(ii) = mean(RD(mask_wm_ct(:)==1));
    rd_st_wm(ii) = mean(RD(mask_wm_st(:)==1));
    ad_wm(ii) = mean(AD(mask_wm_ct(:)==1));
    ad_st_wm(ii) = mean(AD(mask_wm_st(:)==1));
    
    % total kurtosis quantities
    mk_wm(ii) = mean(MKT(mask_wm_ct(:)==1));
    mk_st_wm(ii) = mean(MKT(mask_wm_st(:)==1));
    rk_wm(ii) = mean(RKT(mask_wm_ct(:)==1));
    rk_st_wm(ii) = mean(RKT(mask_wm_st(:)==1));
    ak_wm(ii) = mean(AKT(mask_wm_ct(:)==1));
    ak_st_wm(ii) = mean(AKT(mask_wm_st(:)==1));
    
    % variance kurtosis quantities
    vk_wm(ii) = mean(MKTv(mask_wm_ct(:)==1));
    vk_st_wm(ii) = mean(MKTv(mask_wm_st(:)==1));
    vrk_wm(ii) = mean(RKTv(mask_wm_ct(:)==1));
    vrk_st_wm(ii) = mean(RKTv(mask_wm_st(:)==1));
    vak_wm(ii) = mean(AKTv(mask_wm_ct(:)==1));
    vak_st_wm(ii) = mean(AKTv(mask_wm_st(:)==1));
    
    % microscopic kurtosis quantities
    uk_wm(ii) = mean(KINTRA(mask_wm_ct(:)==1));
    uk_st_wm(ii) = mean(KINTRA(mask_wm_st(:)==1));
    urk_wm(ii) = mean(RKTi(mask_wm_ct(:)==1));
    urk_st_wm(ii) = mean(RKTi(mask_wm_st(:)==1));
    uak_wm(ii) = mean(AKTi(mask_wm_ct(:)==1));
    uak_st_wm(ii) = mean(AKTi(mask_wm_st(:)==1));
    
    % CTI variance quantities kaniso/Kiso/Css
    kaniso_wm(ii) = mean(KANISO(mask_wm_ct(:)==1));
    kaniso_st_wm(ii) = mean(KANISO(mask_wm_st(:)==1));
    kiso_wm(ii) = mean(KISO(mask_wm_ct(:)==1));
    kiso_st_wm(ii) = mean(KISO(mask_wm_st(:)==1));
    cov_wm(ii) = mean(Cov_ss_n(mask_wm_ct(:)==1));
    cov_st_wm(ii) = mean(Cov_ss_n(mask_wm_st(:)==1));
    
    % micro-anisotropic/dispersion quantities
    op_wm(ii) = mean(OP(mask_wm_ct(:)==1));
    op_st_wm(ii) = mean(OP(mask_wm_st(:)==1));
    ufa_wm(ii) =  mean(uFA(mask_wm_ct(:)==1));
    ufa_st_wm(ii) = mean(uFA(mask_wm_st(:)==1));
    
end

%% compute the averages

% diffusion quantities
mfa_wm = mean(fa_wm);
mfa_st_wm = mean(fa_st_wm);
mmd_wm = mean(md_wm);
mmd_st_wm = mean(md_st_wm);
mrd_wm = mean(rd_wm);
mrd_st_wm = mean(rd_st_wm);
mad_wm = mean(ad_wm);
mad_st_wm = mean(ad_st_wm);

% total kurtosis quantities
mmk_wm = mean(mk_wm);
mmk_st_wm = mean(mk_st_wm);
mrk_wm = mean(rk_wm);
mrk_st_wm = mean(rk_st_wm);
mak_wm = mean(ak_wm);
mak_st_wm = mean(ak_st_wm);

% variance kurtosis quantities
mvk_wm = mean(vk_wm);
mvk_st_wm = mean(vk_st_wm);
mvrk_wm = mean(vrk_wm);
mvrk_st_wm = mean(vrk_st_wm);
mvak_wm = mean(vak_wm);
mvak_st_wm = mean(vak_st_wm);

% microscopic kurtosis quantities
muk_wm = mean(uk_wm);
muk_st_wm = mean(uk_st_wm);
murk_wm = mean(urk_wm);
murk_st_wm = mean(urk_st_wm);
muak_wm = mean(uak_wm);
muak_st_wm = mean(uak_st_wm);

% CTI variance quantities kaniso/Kiso/Css
mkaniso_wm = mean(kaniso_wm);
mkaniso_st_wm = mean(kaniso_st_wm);
mkiso_wm = mean(kiso_wm);
mkiso_st_wm = mean(kiso_st_wm);
mcov_wm = mean(cov_wm);
mcov_st_wm = mean(cov_st_wm);

% micro-anisotropic/dispersion quantities
mop_wm  = mean(op_wm);
mop_st_wm  = mean(op_st_wm);
mufa_wm  = mean(ufa_wm);
mufa_st_wm  = mean(ufa_st_wm);

%% Compute relative differences between ct vs st metrics by:
% (Mst-Mct)/Mct = Mst/Mct - 1;

% diffusion quantities
delta_fa = fa_st_wm ./ fa_wm - 1;
delta_md = md_st_wm ./ md_wm - 1;
delta_rd = rd_st_wm ./ rd_wm - 1;
delta_ad = ad_st_wm ./ ad_wm - 1;

% total kurtosis quantities
delta_mk = mk_st_wm ./ mk_wm - 1;
delta_rk = rk_st_wm ./ rk_wm - 1;
delta_ak = ak_st_wm ./ ak_wm - 1;

% variance kurtosis quantities
delta_vk = vk_st_wm ./ vk_wm - 1;
delta_vrk = vrk_st_wm ./ vrk_wm - 1;
delta_vak = vak_st_wm ./ vak_wm - 1;

% microscopic kurtosis quantities
delta_uk = uk_st_wm ./ uk_wm - 1;
delta_urk = urk_st_wm ./ urk_wm - 1;
delta_uak = uak_st_wm ./ uak_wm - 1;

% CTI variance quantities kaniso/Kiso/Css
delta_kaniso = kaniso_st_wm ./ kaniso_wm - 1;
delta_kiso = kiso_st_wm ./ kiso_wm - 1;
delta_cov = cov_st_wm ./ cov_wm - 1;

%% Compute relative differences between the MEAN of ct vs st metrics:

% Diffusion quantities
delta_mfa = mfa_st_wm ./ mfa_wm - 1
delta_mmd = mmd_st_wm ./ mmd_wm - 1
delta_mrd = mrd_st_wm ./ mrd_wm - 1
delta_mad = mad_st_wm ./ mad_wm - 1

% Total kurtosis quantities
delta_mmk = mmk_st_wm ./ mmk_wm - 1;
delta_mrk = mrk_st_wm ./ mrk_wm - 1;
delta_mak = mak_st_wm ./ mak_wm - 1;

% Variance kurtosis quantities
delta_mvk = mvk_st_wm ./ mvk_wm - 1;
delta_mvrk = mvrk_st_wm ./ mvrk_wm - 1;
delta_mvak = mvak_st_wm ./ mvak_wm - 1;

% Microscopic kurtosis quantities
delta_muk = muk_st_wm ./ muk_wm - 1;
delta_murk = murk_st_wm ./ murk_wm - 1;
delta_muak = muak_st_wm ./ muak_wm - 1;

% CTI variance quantities (Kaniso/Kiso/Css)
delta_mkaniso = mkaniso_st_wm ./ mkaniso_wm - 1;
delta_mkiso = mkiso_st_wm ./ mkiso_wm - 1;
delta_mcov = mcov_st_wm ./ mcov_wm - 1;