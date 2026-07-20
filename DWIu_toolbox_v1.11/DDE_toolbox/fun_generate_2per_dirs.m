function [dir1, dir2] = fun_generate_2per_dirs(V)
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

alfa = 0;% half of points because is symetric

[Vx, Vy] = size(V);
if Vx < Vy
    V = V';
    Nv = Vy;
else
    Nv = Vx;
end
per_shift = 2*pi/Nv; % Just to give some more homogenety in
% prependicular directions
%
% If this shift was not done, one direction will
% be always aligned to x=0 plane.
dir1 = zeros(Nv, 3);
dir2 = zeros(Nv, 3);

%dir1(1:Nv, :) = V;
%dir2(1:Nv, :) = V;

for di = 1:Nv
    dir_pa = V(di, :);
    calfa = cos(alfa)';
    salfa = sin(alfa)';
    alfa = alfa + per_shift;
    if(abs(dir_pa(1))>0.999)
        Dir_per = Perp_circle_y(dir_pa(1),dir_pa(2),dir_pa(3),1,calfa,salfa);
    else
        Dir_per = Perp_circle_x(dir_pa(1),dir_pa(2),dir_pa(3),1,calfa,salfa);
    end
    dir1(di, :) = Dir_per;
    dir2(di, :) = cross(dir_pa, Dir_per);
    %disp(alfa)
end

