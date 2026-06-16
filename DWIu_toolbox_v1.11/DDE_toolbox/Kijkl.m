function Kijkl=Kijkl(AllD,DT,f,i,j,k,l)
% Reconstruct Z tensor element (i,j,k,l)
%
% Inputs
% ------
% AllD: matrix (3x3xNcomps)
%    Matrix containing diffusion tensors of microenviroments
% DT: matrix (3x3)
%    Total diffusion tensor
% f: vector (Ncomps)
%    volume fractions
% i: int 
%    Element of tensor 1st dimention
% j: int 
%    Element of tensor 2nd dimention
% k: int 
%    Element of tensor 3rd dimention
% l: int 
%    Element of tensor 4th dimention
%
% Outputs
% -------
% Kijkl: float 
%     Z tensor element (i,j,k,l)
%
% Reference
% ---------
%
%
% Developer
% ---------
% Rafael Neto Henriques 29/06/2018

ncomp = size(AllD,3);
meanDijDkl = 0;

for comp=1:ncomp
    Dc=squeeze(AllD(:,:,comp));
    meanDijDkl = meanDijDkl+f(1,comp)*Dc(i,j)*Dc(k,l);
end
Kijkl = meanDijDkl - DT(i,j)*DT(k,l);
