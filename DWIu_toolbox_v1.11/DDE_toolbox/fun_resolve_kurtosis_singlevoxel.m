function [ktotal, kaniso, kiso, kintra, Convert, MD] = fun_resolve_kurtosis_singlevoxel(pars)

D11 = pars(2);
D22 = pars(3);
D33 = pars(4);
D13 = pars(5);
D12 = pars(6);
D23 = pars(7);

MD = (D11 + D22 + D33) / 3;

W1111 = pars(8);
W2222 = pars(9);
W3333 = pars(10);
W1122 = pars(17);
W1133 = pars(18);
W2233 = pars(19);

C1111 = pars(23);
C2222 = pars(24);
C3333 = pars(25);
C1122 = pars(32);
C1133 = pars(33);
C2233 = pars(34);
C1212 = pars(38);
C1313 = pars(39);
C2323 = pars(40);

% Compute total variance
MW = 1/5 * (W1111 + W2222 + W3333 + 2*W1122 + 2*W1133 + 2*W2233);
MW = MW .* MD .* MD;

Dsumele = D11.^2 + D22.^2 + D33.^2 + 2*D12.^2 + 2*D13.^2 + 2*D23.^2;
Convert = 2/5 * Dsumele - 6/5 * (MD.^2);

TV = MW + Convert;

% Compute variance of tensor magnitudes
VMD = 1/9 * (C1111 + C2222 + C3333 + 2*C1122 + 2*C1133 + 2*C2233);

% Compute averaged eigenvalue variance
EVL = 2/9*(C1111+D11.^2 + C2222+D22.^2 + C3333+D33.^2 - ...
    (C1122+D11.*D22) - (C1133+D11.*D33) - (C2233+D22.*D33) + ...
    3*(C1212+D12.^2 + C1313+D13.^2 + C2323+D23.^2));

% Compute different kurtosis sources
ktotal = TV ./ (MD.^2);
kaniso = 1.2*(EVL) ./ (MD.^2);
kiso = 3*VMD ./ (MD.^2);

ktotal(MD==0) = 0;
kaniso(MD==0) = 0;
kiso(MD==0) = 0;

kintra = ktotal - kaniso - kiso;


