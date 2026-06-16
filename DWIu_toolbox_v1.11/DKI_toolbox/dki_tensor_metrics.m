function [MK, AK, RK] = dki_tensor_metrics(DT, KT, mask)
% Computes standard DKI metrics from diffusion and kurtosis tensors
%
%     Parameters
%     ----------
%     DT : array (..., 6)
%         Array containing the elements of the diffusion tensor in the
%         following order: Dxx, Dyy, Dzz, Dxy, Dxz, Dyz
%     KT : array (..., 15)
%         Array containing the elements of the diffusion tensor in the
%         following order: Wxxxx, Wyyyy, Wzzzz, Wxxxy, Wxxxz, wxyyy
%         Wyyyz, Wxzzz, Wyzzz, Wxxyy, Wxxzz, Wyyzz, Wxxyz, xyzz
%     mask : array (..., )
%         Array containing true values for voxels to be processed
%     kbounds : array (..., ), optional
%         Bounds of kurtosis. Directional kurtosis out of this bounds will
%         be truncated. Default = [-1, 10]
%
%     Returns
%     -------
%     MK :
%         Mean tensor kurtosis
%     AK :
%         Axial tensor kurtosis
%     RK :
%         Radial tensor kurtosis


SIZ = size(DT);
dformat = length(SIZ);

if dformat == 4 % volume
    [MK, AK, RK] = volume_kmeasures(DT, KT, mask, SIZ);
elseif dformat == 3 % slice
    [MK, AK, RK] =  slice_kmeasures(DT, KT, mask, SIZ);
elseif dformat == 2% line scanning
    [MK, AK, RK] = line_kmeasures(DT, KT, mask, SIZ);
else % single voxel
    [MK, AK, RK] = voxel_kmeasures(DT, KT, mask);
end
end

function [mk, ak, rk] = volume_kmeasures(dt, kt, mask, N)
for z = N(3):-1:1
    [mkz, akz, rkz] = slice_kmeasures(squeeze(dt(:, :, z, :)), ...
        squeeze(kt(:, :, z, :)), ...
        squeeze(mask(:, :, z)), N);
    mk(:, :, z) = mkz;
    ak(:, :, z) = akz;
    rk(:, :, z) = rkz;
    disp(['slice N=', num2str(z)])
end
end

function [mk, ak, rk] = slice_kmeasures(dt, kt, mask, N)
for y = N(2):-1:1
    [mky, aky, rky] = line_kmeasures(squeeze(dt(:, y, :)), ...
        squeeze(kt(:, y, :)), squeeze(mask(:, y)), N);
    mk(:, y) = mky;
    ak(:, y) = aky;
    rk(:, y) = rky;
end
end

function [mk, ak, rk] = line_kmeasures(dt, kt, mask, N)
for x = N(1):-1:1
    [mkx, akx, rkx] = voxel_kmeasures(dt(x, :), kt(x, :), ...
        mask(x));
    mk(x) = mkx;
    ak(x) = akx;
    rk(x) = rkx;
    
end
end

function [mk, ak, rk] = voxel_kmeasures(dt, kt, mask)
if mask == 1
    
    Dv=[dt(1) dt(4) dt(5);...
        dt(4) dt(2) dt(6);...
        dt(5) dt(6) dt(3)];
    
    [Vecs,L]=eig(Dv);
    dL=diag(L);
    
    [sL, is]=sort(dL, 'descend');
    sVecs = Vecs(:, is);
    
    md = mean(sL);
    ad = sL(1);
    rd = (sL(2) + sL(3))/2;
    
    W1111 = Wrotate(kt, 1, 1, 1, 1, sVecs);
    W2222 = Wrotate(kt, 2, 2, 2, 2, sVecs);
    W3333 = Wrotate(kt, 3, 3, 3, 3, sVecs);
    W2233 = Wrotate(kt, 2, 2, 3, 3, sVecs);
    
    %mean kurtosis tensor
    mk = (kt(1) + kt(2) + kt(3) + 2*kt(10) + 2*kt(11) + 2*kt(12))/5;
    
    % radial kurtosis tensor
    rk = 3/8 * (W2222 + W3333 + 2*W2233);
    rk = rk * md^2 / (rd^2);
    
    % axial kurtosis tensor
    ak = W1111 * md^2 / (ad^2);

else
    mk = 0;
    ak = 0;
    rk = 0;
end
end

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
end

