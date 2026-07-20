function [D11,D22,D33,D12,D13,D23,...
          W1111, W2222, W3333, W1112, W1113,...
          W1222, W2223, W1333, W2333, W1122,...
          W1133, W2233, W1123, W1223, W1233,...
          K1111, K2222, K3333, K1112, K1113, ...
          K1222, K2223, K1333, K2333, K1122, ...
          K1133, K2233, K1123, K1322, K1233, ...
          K1212, K1313, K2323, K1223, K1323, K1213, S0, SSE]=...
fun_DDE_ULLS_comp_b2_rh(data_in,data_mask, bval1, bval2, bvec1, bvec2)
%% ULLS DDE
% Implemented by Rafael Neto Henriques
% June 2018
%%

[Nx, Ny, Nz, Nvol]=size(data_in);


% Minimize ||AX-B||^2
% Where A, B are:

Ad = zeros(Nvol, 6);
Ak = zeros(Nvol, 15);
Az = zeros(Nvol, 21);

for v=1:Nvol
    b1 = bval1(v);
    b2 = bval2(v);
    Ad(v,1:6)=[b1*bvec1(1,v)^2 + b2*bvec2(1,v)^2, ...
               b1*bvec1(2,v)^2 + b2*bvec2(2,v)^2, ...
               b1*bvec1(3,v)^2 + b2*bvec2(3,v)^2, ..., ...
               2*b1*bvec1(1,v)*bvec1(2,v) + 2*b2*bvec2(1,v)*bvec2(2,v),...
               2*b1*bvec1(1,v)*bvec1(3,v) + 2*b2*bvec2(1,v)*bvec2(3,v), ...
               2*b1*bvec1(2,v)*bvec1(3,v) + 2*b2*bvec2(2,v)*bvec2(3,v)];
    
    Ak(v,1:15)=[b1*b1*bvec1(1,v)^4 + b2*b2*bvec2(1,v)^4,... %xxxx
        b1*b1*bvec1(2,v)^4 + b2*b2*bvec2(2,v)^4, ... %yyyy
        b1*b1*bvec1(3,v)^4 + b2*b2*bvec2(3,v)^4, ... %zzzz
        4*b1*b1*bvec1(1,v)^3*bvec1(2,v) + 4*b2*b2*bvec2(1,v)^3*bvec2(2,v),... %xxxy
        4*b1*b1*bvec1(1,v)^3*bvec1(3,v) + 4*b2*b2*bvec2(1,v)^3*bvec2(3,v),... %xxxz
        4*b1*b1*bvec1(2,v)^3*bvec1(1,v) + 4*b2*b2*bvec2(2,v)^3*bvec2(1,v),... %yyyx
        4*b1*b1*bvec1(2,v)^3*bvec1(3,v) + 4*b2*b2*bvec2(2,v)^3*bvec2(3,v),... %yyyz
        4*b1*b1*bvec1(3,v)^3*bvec1(1,v) + 4*b2*b2*bvec2(3,v)^3*bvec2(1,v),... %zzzx
        4*b1*b1*bvec1(3,v)^3*bvec1(2,v) + 4*b2*b2*bvec2(3,v)^3*bvec2(2,v),... %zzzy
        6*b1*b1*bvec1(1,v)^2*bvec1(2,v)^2 + 6*b2*b2*bvec2(1,v)^2*bvec2(2,v)^2,... %xxyy
        6*b1*b1*bvec1(1,v)^2*bvec1(3,v)^2 + 6*b2*b2*bvec2(1,v)^2*bvec2(3,v)^2,... %xxzz
        6*b1*b1*bvec1(2,v)^2*bvec1(3,v)^2 + 6*b2*b2*bvec2(2,v)^2*bvec2(3,v)^2,... %yyzz
        12*b1*b1*bvec1(1,v)^2*bvec1(2,v)*bvec1(3,v) + 12*b2*b2*bvec2(1,v)^2*bvec2(2,v)*bvec2(3,v),... %xxyz
        12*b1*b1*bvec1(2,v)^2*bvec1(1,v)*bvec1(3,v) + 12*b2*b2*bvec2(2,v)^2*bvec2(1,v)*bvec2(3,v),... %yyxz
        12*b1*b1*bvec1(3,v)^2*bvec1(1,v)*bvec1(2,v) + 12*b2*b2*bvec2(3,v)^2*bvec2(1,v)*bvec2(2,v)]; %zzxy

    Az(v, 1:21) = [b1*bvec1(1,v)^2 * b2*bvec2(1,v)^2, ... %xxxx
                   b1*bvec1(2,v)^2 * b2*bvec2(2,v)^2, ... %yyyy
                   b1*bvec1(3,v)^2 * b2*bvec2(3,v)^2, ... %zzzz
                   2 * (b1*bvec1(1,v)^2 * b2*bvec2(1,v)*bvec2(2,v) + ...
                        b1*bvec1(1,v)*bvec1(2,v) * b2*bvec2(1,v)^2), ... %xxxy
                   2 * (b1*bvec1(1,v)^2 * b2*bvec2(1,v)*bvec2(3,v) + ...
                        b1*bvec1(1,v)*bvec1(3,v) * b2*bvec2(1,v)^2), ... %xxxz
                   2 * (b1*bvec1(2,v)^2 * b2*bvec2(2,v)*bvec2(1,v) + ...
                        b1*bvec1(2,v)*bvec1(1,v) * b2*bvec2(2,v)^2), ... %xyyy
                   2 * (b1*bvec1(2,v)^2 * b2*bvec2(2,v)*bvec2(3,v) + ...
                        b1*bvec1(2,v)*bvec1(3,v) * b2*bvec2(2,v)^2), ... %yyyz (7)
                   2 * (b1*bvec1(3,v)^2 * b2*bvec2(3,v)*bvec2(1,v) + ...
                        b1*bvec1(3,v)*bvec1(1,v) * b2*bvec2(3,v)^2), ... %xzzz (8)
                   2 * (b1*bvec1(3,v)^2 * b2*bvec2(3,v)*bvec2(2,v) + ...
                        b1*bvec1(3,v)*bvec1(2,v) * b2*bvec2(3,v)^2), ... %zzzy (9)
                   b1*bvec1(1,v)^2 * b2*bvec2(2,v)^2 + ...
                       b1*bvec1(2,v)^2 * b2*bvec2(1,v)^2 , ... %xxyy (10)
                   b1*bvec1(1,v)^2 * b2*bvec2(3,v)^2 + ...
                       b1*bvec1(3,v)^2 * b2*bvec2(1,v)^2 , ... %xxzz (11)
                   b1*bvec1(2,v)^2 * b2*bvec2(3,v)^2 + ...
                       b1*bvec1(3,v)^2 * b2*bvec2(2,v)^2 , ... %yyzz (12)
                   2 * (b1*bvec1(1,v)^2 * b2*bvec2(2,v)*bvec2(3,v) + ...
                        b1*bvec1(2,v)*bvec1(3,v) * b2*bvec2(1,v)^2), ... %xxyz (13)
                   2 * (b1*bvec1(2,v)^2 * b2*bvec2(1,v)*bvec2(3,v) + ...
                        b1*bvec1(1,v)*bvec1(3,v) * b2*bvec2(2,v)^2), ... %yyxz (14)
                   2 * (b1*bvec1(3,v)^2 * b2*bvec2(1,v)*bvec2(2,v) + ...
                        b1*bvec1(1,v)*bvec1(2,v) * b2*bvec2(3,v)^2), ... %zzxy (15)
                   4 * (b1*bvec1(1,v)*bvec1(2,v) * b2*bvec2(1,v)*bvec2(2,v)),... %xyxy
                   4 * (b1*bvec1(1,v)*bvec1(3,v) * b2*bvec2(1,v)*bvec2(3,v)),... %xzxz
                   4 * (b1*bvec1(2,v)*bvec1(3,v) * b2*bvec2(2,v)*bvec2(3,v)),... %yzyz
                   4 * (b1*bvec1(1,v)*bvec1(2,v) * b2*bvec2(2,v)*bvec2(3,v) + ...
                        b1*bvec1(2,v)*bvec1(3,v) * b2*bvec2(1,v)*bvec2(2,v)) ... %xyyz
                   4 * (b1*bvec1(1,v)*bvec1(3,v) * b2*bvec2(2,v)*bvec2(3,v) + ...
                        b1*bvec1(2,v)*bvec1(3,v) * b2*bvec2(1,v)*bvec2(3,v)) ... %xzzy
                   4 * (b1*bvec1(1,v)*bvec1(2,v) * b2*bvec2(1,v)*bvec2(3,v) + ...
                        b1*bvec1(1,v)*bvec1(3,v) * b2*bvec2(1,v)*bvec2(2,v)) ... %yxxz
];
end

%A
A=[ -Ad 1/6*Ak Az ones(Nvol, 1)]; %SO is now the last column

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

SSE = Z0;

S0=Z0;

for k=1:Nz
    for j=1:Ny
        for i=1:Nx  
            if(data_mask(i,j,k)==1)  
                %B
                B=log(squeeze(data_in(i,j,k,:)));
                
                indinf=isinf(B);
                if sum(indinf)~=0
                    %minnotinf=min(B(~indinf)); 
                    %maximum diffusion posible take it other 
                    %direction but perhaps other aproaches 
                    %can be better
                    B(indinf)=-10000;  
                end
                
                % ULLS
                piA=pinv(A); %piA pseudoinverse of A
                X=piA*B;
                
                %% Diffusion parameters
                % diffusion tensor
                D11(i,j,k)=X(1);
                D22(i,j,k)=X(2);
                D33(i,j,k)=X(3);
                D12(i,j,k)=X(4);
                D13(i,j,k)=X(5);
                D23(i,j,k)=X(6);
                
                if isnan(X(1))
                    disp('NaN');
                end
                
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
