function [ZT, KT, DT] = Zgenerator_directions(DAi, DRi, DAe, DRe, fie, fp, V)

%addpath(['..', filesep , 'DKI_toolbox'])

Di=[DAi 0 0;0 DRi 0;0 0 DRi];
De=[DAe 0 0;0 DRe 0;0 0 DRe];

NoF = length(fp);

AllDi=zeros(3,3,NoF);
AllDe=AllDi;
DTci=AllDi;
DTce=AllDi;
f_p1=fp; % compartment fractions part 1 (intracelular compartments)
f_p2=fp; % compartment fractions part 2 (extracelular compartments)


for c=1:NoF

        f_p1(c)=fp(c)*fie;
        f_p2(c)=fp(c)*(1-fie);


    Xa=V(c, :);
    
    if Xa(1) == 1
        Xb=cross(Xa,[0 1 0]);
        Xb=Xb/norm(Xb);
        Xc=cross(Xa,Xb);
        Xc=Xc/norm(Xc);
    else
        Xb=cross(Xa,[1 0 0]);
        Xb=Xb/norm(Xb);
        Xc=cross(Xa,Xb);
        Xc=Xc/norm(Xc);
    end
    
    Evector=[Xa',Xb',Xc'];
    
    Dci=Evector*Di*pinv(Evector);
    AllDi(:,:,c)=Dci;
    DTci(:,:,c)=f_p1(1,c)*Dci; % Matrix to later compute DT;
    
    Dce=Evector*De*pinv(Evector);
    AllDe(:,:,c)=Dce;
    DTce(:,:,c)=f_p2(1,c)*Dce; % Matrix to later compute DT;
end

% Compute DT, f and AllD
f=[f_p1 f_p2];
DTc=zeros(3,3,NoF*2);
AllD=DTc;
DTc(:,:,1:NoF)=DTci;
DTc(:,:,NoF+1:2*NoF)=DTce;
DT=sum(DTc,3);
AllD(:,:,1:NoF)=AllDi;
AllD(:,:,NoF+1:2*NoF)=AllDe;

% matrix para plot posterior
Dv=[DT(1,1) DT(2,2) DT(3,3) DT(1,2) DT(1,3) DT(2,3)];

% Calculo do tensor de kurtosis
W = sim_kt(AllD,DT,f);
%WT=Wcons(W);

ZT = sim_zt(AllD,DT,f);

%Ziiii = Kijkl(AllD, DT, f, 1, 1, 1, 1);
%Zjjjj = Kijkl(AllD, DT, f, 2, 2, 2, 2);
%Zzzzz = Kijkl(AllD, DT, f, 3, 3, 3, 3);
%Zxxyy = Kijkl(AllD, DT, f, 1, 1, 2, 2);
%Zxxzz = Kijkl(AllD, DT, f, 1, 1, 3, 3);
%Zyyzz = Kijkl(AllD, DT, f, 2, 2, 3, 3); 
% 1 2 3 10 11 12
%VDiso = (Ziiii + Zjjjj + Zzzzz + 2*Zxxyy + 2*Zxxzz + 2*Zyyzz)/9;

VDiso_gt = 0;
MD = trace(DT)/3;
for c=1:(NoF*2)
    Dc = AllD(:, :, c);
    MDc = trace(Dc)/3;
    VDiso_gt = VDiso_gt + f(c)*(MDc-MD)^2;
end 

DT = Dv;
KT = W;
