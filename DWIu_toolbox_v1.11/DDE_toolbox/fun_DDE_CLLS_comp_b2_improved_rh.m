function [D11,D22,D33,D12,D13,D23,...
    W1111, W2222, W3333, W1112, W1113,...
    W1222, W2223, W1333, W2333, W1122,...
    W1133, W2233, W1123, W1223, W1233,...
    K1111, K2222, K3333, K1112, K1113, ...
    K1222, K2223, K1333, K2333, K1122, ...
    K1133, K2233, K1123, K1322, K1233, ...
    K1212, K1313, K2323, K1223, K1323, K1213, S0, SSE, Ndiff, Nkurt, Ncov]=...
    fun_DDE_CLLS_comp_b2_improved_rh(data_in,data_mask, bval1, bval2, bvec1, bvec2, V)
%% CLLS DDE
% Implemented by Rafael Neto Henriques
% June 2024
%%

[Nx, Ny, Nz, Nvol] = size(data_in);

bt = bval1 + bval2;

% CTI design_matrix
A = fun_CTI_desing_matrix(bval1, bval2, bvec1, bvec2);


%% Contraint design matrix
ndir = size(V, 1);
[V2, V3] = fun_generate_2per_dirs(V);
bdc1 = [V; V2];
bdc2 = [V; V3];
bvc1 = ones(ndir*2, 1)/2;
bvc2 = ones(ndir*2, 1)/2;

% QTI design matrix to contraint covariant tensors
bxx = bvc1 .* bdc1(:, 1).^2 + bvc2 .* bdc2(:, 1).^2;
byy = bvc1 .* bdc1(:, 2).^2 + bvc2 .* bdc2(:, 2).^2;
bzz = bvc1 .* bdc1(:, 3).^2 + bvc2 .* bdc2(:, 3).^2;
bxy = bvc1 .* bdc1(:, 1).* bdc1(:, 2) +...
    bvc2 .* bdc2(:, 1).* bdc2(:, 2);
bxz = bvc1 .* bdc1(:, 1).* bdc1(:, 3) +...
    bvc2 .* bdc2(:, 1).* bdc2(:, 3);
byz = bvc1 .* bdc1(:, 2).* bdc1(:, 3) +...
    bvc2 .* bdc2(:, 2).* bdc2(:, 3);
Aqti = fun_QTI_desing_matrix(bxx, byy, bzz, bxy, bxz, byz);

bdc1 = [V; bdc1];
bdc2 = [V; bdc2];
bvc1 = [ones(ndir, 1); bvc1];
bvc2 = [zeros(ndir, 1); bvc2];
Acti = fun_CTI_desing_matrix(bvc1, bvc2, bdc1, bdc2);

% constraints Kim = 0 (so that I don't need estimation of app_diff apriori

% contrain matrix (diffusivities)
ad = Acti;
ad(:, 7:end) = 0;
%app_diff = -ad * X0;
%app_diff(app_diff<0) = 0;

% contrain matrix (kurtosis)
zzz = zeros(ndir, 1);
%d = Kim .* app_diff.^2;
d = [zzz(:); zzz(:); zzz(:)];
ak = Acti;
ak(:, 1:6) = 0;
ak(:, end) = 0;

% contrain matrix (covariance)
%av = ad(1:2*ndir, :) * 0;
%av(:, 22:42) = Aqti(:, 7:27);

% contrain only trace of covariance tensor
Aiso = fun_QTI_desing_matrix(1/3, 1/3, 1/3, 0, 0, 0);
av_iso = ad(1, :) * 0;
av_iso(22:42) = Aiso(:, 7:27);

% concatenate
Const = [ad; -ak; -av_iso];% -av_iso];
ddd = [zzz(:); zzz(:); zzz(:);...
    -d(:); ...
    0]; %zzz(:); zzz(:)];% 0];

%% 

Z0=zeros(Nx,Ny,Nz);

D11=Z0; D22=Z0; D33=Z0; D12=Z0; D13=Z0; D23=Z0;

W1111=Z0; W2222=Z0; W3333=Z0; W1112=Z0; W1113=Z0;
W1222=Z0; W2223=Z0; W1333=Z0; W2333=Z0; W1122=Z0;
W1133=Z0; W2233=Z0; W1123=Z0; W1223=Z0; W1233=Z0;

K1111=Z0; K2222=Z0; K3333=Z0; K1112=Z0; K1113=Z0;
K1222=Z0; K2223=Z0; K1333=Z0; K2333=Z0; K1122=Z0;
K1133=Z0; K2233=Z0; K1123=Z0; K1322=Z0; K1233=Z0;
K1212=Z0; K1313=Z0; K2323=Z0; K1223=Z0; K1323=Z0;
K1213=Z0;

Ndiff = Z0 + 3*ndir;
Nkurt = Z0 + 3*ndir;
Ncov = Z0 + 1;

SSE = Z0;

S0=Z0;

options = optimoptions('lsqlin','Algorithm','interior-point',...
    'MaxIterations',20,'linearsolver', 'Dense', 'Display','none',...
    'ConstraintTolerance', 1.0000e-08, ...
    'OptimalityTolerance', 1.0000e-08, ...
    'StepTolerance', 1.0000e-12);

for k=1:Nz
    for j=1:Ny
        for i=1:Nx
            if(data_mask(i,j,k)==1)
                %B
                B=log(squeeze(data_in(i,j,k,:)));
                
                indinf=isinf(B);
                if sum(indinf)~=0
                    B(indinf)=-10000;
                end
                
                % ULLS
                piA=pinv(A); %piA pseudoinverse of A
                X=piA*B;
                
                %
                cond1_past = sum(ad * X<=0);
                Ndiff(i,j,k) = cond1_past;
                
                app_kurt = ak * X;
                cond2_past = sum(-app_kurt <=0);
                Nkurt(i,j,k) = cond2_past;
                
                app_cov = av_iso * X;
                cond3_past = sum(-app_cov <=0);
                Ncov(i,j,k) = cond3_past;
                
                if cond1_past + cond2_past + cond3_past < 6*ndir+1
                    X = lsqlin(A, B, Const, ddd, [], [], [], [], [], options);
                end
                
                %% Diffusion parameters
                % diffusion tensor
                D11(i,j,k)=X(1);
                D22(i,j,k)=X(2);
                D33(i,j,k)=X(3);
                D12(i,j,k)=X(4);
                D13(i,j,k)=X(5);
                D23(i,j,k)=X(6);
                
                %if isnan(X(1))
                %    disp('NaN');
                %end
                
                %mean diffusion
                MD=(D11(i,j,k)+D22(i,j,k)+D33(i,j,k))/3.0;
                
                %% t2 imaging
                S0(i,j,k)=exp(X(43));
                
                %% Kurtosis parameters
                
                W=X(7:21)/(MD^2);
                
                W1111(i,j,k)=W(1);
                W2222(i,j,k)=W(2);
                W3333(i,j,k)=W(3);
                W1112(i,j,k)=W(4);
                W1113(i,j,k)=W(5);
                W1222(i,j,k)=W(6);
                W2223(i,j,k)=W(7);
                W1333(i,j,k)=W(8);
                W2333(i,j,k)=W(9);
                W1122(i,j,k)=W(10);
                W1133(i,j,k)=W(11);
                W2233(i,j,k)=W(12);
                W1123(i,j,k)=W(13);
                W1223(i,j,k)=W(14);
                W1233(i,j,k)=W(15);
                
                %% Ztensor parameters
                K1111(i,j,k)=X(22); %1
                K2222(i,j,k)=X(23); %2
                K3333(i,j,k)=X(24); %3
                K1112(i,j,k)=X(25); %4
                K1113(i,j,k)=X(26); %5
                K1222(i,j,k)=X(27); %6
                K2223(i,j,k)=X(28); %7
                K1333(i,j,k)=X(29); %8
                K2333(i,j,k)=X(30); %9
                K1122(i,j,k)=X(31); %10
                K1133(i,j,k)=X(32); %11
                K2233(i,j,k)=X(33); %12
                K1123(i,j,k)=X(34); %13
                K1322(i,j,k)=X(35); %14
                K1233(i,j,k)=X(36); %15
                K1212(i,j,k)=X(37); %16
                K1313(i,j,k)=X(38); %17
                K2323(i,j,k)=X(39); %18
                K1223(i,j,k)=X(40); %19
                K1323(i,j,k)=X(41); %20
                K1213(i,j,k)=X(42); %21
                
                SSE(i,j,k) = sum((squeeze(data_in(i,j,k,:))-exp(A*X)).^2);
                
            end
        end
    end
    disp(k)
end
