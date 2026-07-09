function [d, k, w] = diffusion_in_cylinder(t, R, N, d0)
% EXACT ANALYTICAL SOLUTION FOR RADIAL DIFFUSION IN A CYLINDER
% t  : diffusion time
% R  : cylinder radius
% N  : number of Bessel roots to include (max 20)
% d0 : intrinsic diffusivity inside the cylinder

% Roots of Bessel functions (precomputed for N up to 20)
% a0: First 20 roots of J1(x) = 0
a0 = [ ...
    3.8317059702,  7.0155866698, 10.1734681351, 13.3236919363, ...
    16.4706300508, 19.6158585105, 22.7600843806, 25.9036720876, ...
    29.0468285349, 32.1896799110, 35.3323075501, 38.4747662348, ...
    41.6170942128, 44.7593189982, 47.9014608872, 51.0435351836, ...
    54.1855536411, 57.3275254381, 60.4694578453, 63.6113566985];

% a1: First 20 roots of J1'(x) = 0
a1 = [ ...
    1.8411837813,  5.3314427735,  8.5363163663, 11.7060049026, ...
    14.8635886339, 18.0155278544, 21.1643697405, 24.3113268069, ...
    27.4570514102, 30.6019234891, 33.7460982566, 36.8896944358, ...
    40.0328135889, 43.1755255474, 46.3178877196, 49.4599507914, ...
    52.6017551066, 55.7433346513, 58.8847188725, 62.0259339598];

% a2: First 20 roots of J2'(x) = 0
a2 = [ ...
    3.0542369282,  6.7061331942,  9.9694678231, 13.1703708560, ...
    16.3475223157, 19.5129130095, 22.6681656515, 25.8145558066, ...
    28.9535090407, 32.0860363980, 35.2128795551, 38.3346082260, ...
    41.4516900898, 44.5645391167, 47.6735235885, 50.7789642646, ...
    53.8811440026, 56.9803157297, 60.0767073289, 63.1705193910];

% Note if you want more roots. Use the code below:

%N_max = 20; 
%a0_dynamic = zeros(1, N_max);
%for n = 1:N_max
    % The roots of J_1(x) asymptotically approach (n + 0.25)*pi
%    guess = (n + 0.25) * pi; 
%    a0_dynamic(n) = fzero(@(x) besselj(1, x), guess);
%end


% Restrict N to avoid indexing errors
N = min(N, 20);
a0 = a0(1:N);
a1 = a1(1:N);
a2 = a2(1:N);

% Exponential decay factors for each mode
e0 = exp(-a0.^2 * d0 * t / R^2);
e1 = exp(-a1.^2 * d0 * t / R^2);
e2 = exp(-a2.^2 * d0 * t / R^2);

% Mean squared displacement (1D projection of 2D radial displacement: <x^2>)
% The asymptotic limit as t -> inf is R^2 / 2
C1 = 4 ./ (a1.^2 .* (a1.^2 - 1));
x2 = (R^2 / 2) - R^2 * sum(e1 .* C1);

% Mean quartic displacement (1D projection: <x^4>)
% The asymptotic limit as t -> inf is 5/8 * R^4
B0 = 24 ./ (a0.^4);
B1 = -12 * (3 * a1.^2 - 8) ./ (a1.^4 .* (a1.^2 - 1));
B2 = 12 ./ (a2.^2 .* (a2.^2 - 4));

x4 = (5/8 * R^4) + R^4 * sum(e0 .* B0) + R^4 * sum(e1 .* B1) + R^4 * sum(e2 .* B2);

% Kurtosis and apparent diffusivity calculation
k = x4 / (x2^2) - 3;
d = x2 / (2 * t);

w = d^2 * k;
end