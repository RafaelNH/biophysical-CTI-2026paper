function b = bingham(u, f0, k1, k2, u1, u2)
% Samples the bingham distribution with parameters f0, k1, k2, u1, and u2
% on points u.
%
% Inputs
% ------
% u: matrix N x 3 
%     Matix containing the 3 coordinates of N points
% f0, k1, k2: scalars
% u1, u2: vectors perpendicular to bingham peak

K = [k1; k2];

uu1 = u(:, 1)*u1(1) + u(:, 2)*u1(2) + u(:, 3)*u1(3);
uu2 = u(:, 1)*u2(1) + u(:, 2)*u2(2) + u(:, 3)*u2(3);

b = f0 * exp(- [uu1.^2, uu2.^2] * K);