function B = Bingham_2D_for_ploting(f0, k1, k2, u1, u2, theta, phi)
% Samples the bingham distribution with parameters f0, k1, k2, u1, and u2
% on points u.
%
% Inputs
% ------
% u: matrix N x 3 
%     Matix containing the 3 coordinates of N points
% f0, k1, k2: scalars
% u1, u2: vectors perpendicular to bingham peak
% theta, phi (messgrid variable containing the spherical coordinates of the
% directions you want to sample)

[x, y, z] = sph2cart(theta, phi, 1);
N = size(x);
na = N(1);
nb = N(2);

B = zeros(N);

for a=1:na
    for b=1:nb
        gnz=[x(a,b), y(a,b), z(a,b)]; 
        B(a, b) = bingham(gnz, f0, k1, k2, u1, u2);
    end
end