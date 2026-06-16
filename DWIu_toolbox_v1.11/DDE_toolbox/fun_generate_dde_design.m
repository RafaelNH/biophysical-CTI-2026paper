function [dir1, dir2] = fun_generate_dde_design(V, nodir, angle)
% Generates a DDE desing (dir1 and dir2) given the parallel directions.
%
% Inputs:
% -------
% V: Matrix (3, Nv)
%    A matrix containing the parallel directions
% nodir: int
%    Number of perpendicular directions per parallel directions
% angle: float (rad)
%    Angle between dir1 and dir2
%
% Output:
% -------
% dir1: Matrix (3, Nv * (nodir+1))
%     A matrix containing DDE first directions
% dir2: Matrix (3, Nv * (nodir+1))
%     A matrix containing DDE second directions
%
% Implemented by
% --------------
% Rafael Neto Henriques (16/04/2018)
%
% Notes
% -----
% Directions are generated in the following order:
% Parallel d 1, 
% Perpendicular d1 1,
% Perpendicular d1 2,
% (...)
% Perpendicular d1 nodir
% Parallel d 2
% Perpendicular d2 1,
% (...)
% Perpendicular d2 nodir
% (...)
% Parallel d(length(V)) 
% (...)
% Perpendicular d(length(V)) nodir

alfa = 0:2*pi/nodir:(2*pi-2*pi/nodir);% half of points because is symetric

[Vx, Vy] = size(V);
if Vx < Vy
    V = V';
    Nv = Vy;
else
    Nv = Vx;
end
per_shift = 2*pi/(nodir*Nv); % Just to give some more homogenety in 
                           % prependicular directions
                           %
                           % If this shift was not done, one direction will
                           % be always aligned to x=0 plane.
dir1 = zeros(Nv * (nodir+1), 3);
dir2 = zeros(Nv * (nodir+1), 3);

%dir1(1:Nv, :) = V;
%dir2(1:Nv, :) = V;

a = 0;
for di = 1:Nv
    dir_pa = V(di, :);        
    calfa = cos(alfa)';
    salfa = sin(alfa)';
    alfa = alfa + per_shift;
    if(abs(dir_pa(1))>0.999)
        Dir_per = Perp_circle_y(dir_pa(1),dir_pa(2),dir_pa(3),nodir,calfa,salfa);
    else
        Dir_per = Perp_circle_x(dir_pa(1),dir_pa(2),dir_pa(3),nodir,calfa,salfa);
    end
    a = a + 1;
    dir1(a, :) = dir_pa;
    dir2(a, :) = dir_pa;
    for r = 1:nodir
        a = a + 1;
        dir1(a, :) = dir_pa;
        dir2(a, :) = Dir_per(r, :);
    end
    disp(alfa)
end

dir3 = dir1 * cos(angle) + dir2 * sin(angle);
nd = sqrt(diag(dir3*dir3'));
dir2 = [dir3(:, 1)./nd, dir3(:, 2)./nd, dir3(:, 3)./nd];