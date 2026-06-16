function [MD, AD, RD, FA, V1, V2, V3, L2, L3] = dti_metrics(DT, mask)
% Computes standard metrics of the diffusion tensor
%
%     Parameters
%     ----------
%     DT : array (..., 6)
%         Array containing the elements of the diffusion tensor in the
%         following order: Dxx, Dyy, Dzz, Dxy, Dxz, Dyz
%     mask : array (..., )
%         Array containing true values for voxels to be processed
%
%     Returns
%     -------
%     MD : array (...)
%         Mean diffusivity
%     AD : array (...)
%         Axial diffusivity
%     RD : array (...)
%         Radial diffusivity
%     FA : array (...)
%         Fractional anisotropy
%     V1 : array (..., 3)
%         Diffusion tensor principal axis
%     V2 : array (..., 3)
%         Diffusion secondary axis
%     V3 : array (..., 3)
%         Diffusion axis of diffusion tensor lower diffusivuty
%     L2 : array (...)
%         Diffusivity along secondary diffusion tensor secondary axis
%     L3 : array (...)
%         Diffusion tensor lower diffusivuty

SIZ = size(DT);
dformat = length(SIZ);

if dformat == 4 % volume
    [MD, AD, RD, FA, V1, V2, V3, L2, L3] = volume_dmeasures(DT, mask, SIZ);
elseif dformat == 3 % slice
    [MD, AD, RD, FA, V1, V2, V3, L2, L3] =  slice_dmeasures(DT, mask, SIZ);
elseif dformat == 2% line scanning
    [MD, AD, RD, FA, V1, V2, V3, L2, L3] = line_dmeasures(DT, mask, SIZ);
else % single voxel
    [MD, AD, RD, FA, V1, V2, V3, L2, L3] = voxel_dmeasures(DT, mask);
end
end

function [md, ad, rd, fa, v1, v2, v3, l2, l3] = volume_dmeasures(dt, mask, N)
for z = N(3):-1:1;
    [mdz, adz, rdz, faz, v1z, v2z, v3z, l2z, l3z] = slice_dmeasures(squeeze(dt(:, :, z, :)), ...
        squeeze(mask(:, :, z)), N);
    md(:, :, z) = mdz;
    ad(:, :, z) = adz;
    rd(:, :, z) = rdz;
    fa(:, :, z) = faz;
    v1(:, :, z, :) = v1z;
    v2(:, :, z, :) = v2z;
    v3(:, :, z, :) = v3z;
    l2(:, :, z) = l2z;
    l3(:, :, z) = l3z;
    disp(['slice N=', num2str(z)])
end
end

function [md, ad, rd, fa, v1, v2, v3, l2, l3] = slice_dmeasures(dt, mask, N)
for y = N(2):-1:1;
    [mdy, ady, rdy, fay, v1y, v2y, v3y, l2y, l3y] = line_dmeasures(squeeze(dt(:, y, :)), ...
        squeeze(mask(:, y)), N);
    md(:, y) = mdy;
    ad(:, y) = ady;
    rd(:, y) = rdy;
    fa(:, y) = fay;
    v1(:, y, :) = v1y;
    v2(:, y, :) = v2y;
    v3(:, y, :) = v3y;
    l2(:, y) = l2y;
    l3(:, y) = l3y;
end
end

function [md, ad, rd, fa, v1, v2, v3, l2, l3] = line_dmeasures(dt, mask, N)
for x = N(1):-1:1;
    [mdx, adx, rdx, fax, v1x, v2x, v3x, l2x, l3x] = voxel_dmeasures(dt(x, :), mask(x));
    md(x) = mdx;
    ad(x) = adx;
    rd(x) = rdx;
    fa(x) = fax;
    v1(x, :) = v1x;
    v2(x, :) = v2x;
    v3(x, :) = v3x;
    l2(x) = l2x;
    l3(x) = l3x;
end
end

function [md, ad, rd, fa, v1, v2, v3, l2, l3] = voxel_dmeasures(dt, mask)
if mask == 1
    % Tensor
    Dv=[dt(1) dt(4) dt(5);...
        dt(4) dt(2) dt(6);...
        dt(5) dt(6) dt(3)];
    
    % Eigen-decomposition
    [Vecs, L]=eig(Dv);
    dL=diag(L);
    [sL, is]=sort(dL);
    l1=sL(3);
    l2=sL(2);
    l3=sL(1);
    
    % DTI measures
    
    md=(l1+l2+l3)/3;
    ad = l1;
    rd = (l2 + l3)/2 ;
    fa=sqrt(3/2*((l1-md)^2+(l2-md)^2+(l3-md)^2 )/(l1^2+l2^2+l3^2) );
    v1=Vecs(:, is(3));
    v2=Vecs(:, is(2));
    v3=Vecs(:, is(1));
    
else
    md = 0;
    ad = 0;
    rd = 0;
    fa = 0;
    v1 = [0, 0, 0];
    v2 = [0, 0, 0];
    v3 = [0, 0, 0];
    l2 = 0;
    l3 = 0;
end
end


