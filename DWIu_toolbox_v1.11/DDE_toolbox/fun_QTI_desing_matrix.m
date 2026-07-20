function A = fun_QTI_desing_matrix(bxx, byy, bzz, bxy, bxz, byz)
% Generate QTI design matrix
%
% Implemented by Rafael Neto Henriques
% May 2024
%%


bxx = bxx(:);
byy = byy(:);
bzz = bzz(:);
bxy = bxy(:);
bxz = bxz(:);
byz = byz(:);

Ad = [bxx, byy, bzz, 2*bxy, 2*bxz, 2*byz];

Az = [bxx.*bxx, byy.*byy, bzz.*bzz, ...
      4*bxx.*bxy, 4*bxx.*bxz, 4*bxy.*byy, 4*byy.*byz, 4*bxz.*bzz, 4*bzz.*byz, ...
      2*bxx.*byy, 2*bxx.*bzz, 2*byy.*bzz,...
      4*bxx.*byz, 4*byy.*bxz, 4*bzz.*bxy, 4*bxy.*bxy, 4*bxz.*bxz, 4*byz.*byz, ...
      8*bxy.*byz, 8*bxz.*byz, 8*bxy.*bxz];

A=[-Ad Az/2 ones(length(bxx), 1)];





