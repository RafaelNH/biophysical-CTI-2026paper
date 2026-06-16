function Wrotated = Wrotate(W, indi, indj, indk, indl, R)
% Rotate kurtosis elements to frame of reference R
% assumes that W = W1111 W2222 W3333 W1112 W1113
%                  W1222 W2223 W1333 W2333 W1122
%                  W1133 W2233 W1123 W1223 W1233
% Rafael Neto Henriques

xyz = 1:3;
Wrotated = 0;

comb=[1 16 81 2 3 8 24 27 54 4 9 36 6 12 18];

for ii = xyz
    for jj = xyz
        for kk = xyz
            for ll = xyz
                B = R(ii, indi)*R(jj, indj)*R(kk, indk)*R(ll, indl);
                Wrotated = Wrotated + B*W(comb==ii*jj*kk*ll);
            end
        end
    end
end


