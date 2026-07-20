function A = fun_CTI_desing_matrix(bval1, bval2, bvec1, bvec2)
% Generate CTI design matrix
%
% Implemented by Rafael Neto Henriques
% May 2024
%%

siz = size(bvec1);
if siz(1)>siz(2)
    bvec1 = bvec1';
    bvec2 = bvec2';
end

Nvol = length(bval1);

Ad = zeros(Nvol, 6);
Ak = zeros(Nvol, 15);
Az = zeros(Nvol, 21);

for v=1:Nvol
    b1 = bval1(v);
    b2 = bval2(v);
    Ad(v,1:6)=[b1*bvec1(1,v)^2 + b2*bvec2(1,v)^2, ...
               b1*bvec1(2,v)^2 + b2*bvec2(2,v)^2, ...
               b1*bvec1(3,v)^2 + b2*bvec2(3,v)^2, ..., ...
               2*b1*bvec1(1,v)*bvec1(2,v) + 2*b2*bvec2(1,v)*bvec2(2,v),...
               2*b1*bvec1(1,v)*bvec1(3,v) + 2*b2*bvec2(1,v)*bvec2(3,v), ...
               2*b1*bvec1(2,v)*bvec1(3,v) + 2*b2*bvec2(2,v)*bvec2(3,v)];
    
    Ak(v,1:15)=[b1*b1*bvec1(1,v)^4 + b2*b2*bvec2(1,v)^4,... %xxxx
        b1*b1*bvec1(2,v)^4 + b2*b2*bvec2(2,v)^4, ... %yyyy
        b1*b1*bvec1(3,v)^4 + b2*b2*bvec2(3,v)^4, ... %zzzz
        4*b1*b1*bvec1(1,v)^3*bvec1(2,v) + 4*b2*b2*bvec2(1,v)^3*bvec2(2,v),... %xxxy
        4*b1*b1*bvec1(1,v)^3*bvec1(3,v) + 4*b2*b2*bvec2(1,v)^3*bvec2(3,v),... %xxxz
        4*b1*b1*bvec1(2,v)^3*bvec1(1,v) + 4*b2*b2*bvec2(2,v)^3*bvec2(1,v),... %yyyx
        4*b1*b1*bvec1(2,v)^3*bvec1(3,v) + 4*b2*b2*bvec2(2,v)^3*bvec2(3,v),... %yyyz
        4*b1*b1*bvec1(3,v)^3*bvec1(1,v) + 4*b2*b2*bvec2(3,v)^3*bvec2(1,v),... %zzzx
        4*b1*b1*bvec1(3,v)^3*bvec1(2,v) + 4*b2*b2*bvec2(3,v)^3*bvec2(2,v),... %zzzy
        6*b1*b1*bvec1(1,v)^2*bvec1(2,v)^2 + 6*b2*b2*bvec2(1,v)^2*bvec2(2,v)^2,... %xxyy
        6*b1*b1*bvec1(1,v)^2*bvec1(3,v)^2 + 6*b2*b2*bvec2(1,v)^2*bvec2(3,v)^2,... %xxzz
        6*b1*b1*bvec1(2,v)^2*bvec1(3,v)^2 + 6*b2*b2*bvec2(2,v)^2*bvec2(3,v)^2,... %yyzz
        12*b1*b1*bvec1(1,v)^2*bvec1(2,v)*bvec1(3,v) + 12*b2*b2*bvec2(1,v)^2*bvec2(2,v)*bvec2(3,v),... %xxyz
        12*b1*b1*bvec1(2,v)^2*bvec1(1,v)*bvec1(3,v) + 12*b2*b2*bvec2(2,v)^2*bvec2(1,v)*bvec2(3,v),... %yyxz
        12*b1*b1*bvec1(3,v)^2*bvec1(1,v)*bvec1(2,v) + 12*b2*b2*bvec2(3,v)^2*bvec2(1,v)*bvec2(2,v)]; %zzxy

    Az(v, 1:21) = [b1*bvec1(1,v)^2 * b2*bvec2(1,v)^2, ... %xxxx
                   b1*bvec1(2,v)^2 * b2*bvec2(2,v)^2, ... %yyyy
                   b1*bvec1(3,v)^2 * b2*bvec2(3,v)^2, ... %zzzz
                   2 * (b1*bvec1(1,v)^2 * b2*bvec2(1,v)*bvec2(2,v) + ...
                        b1*bvec1(1,v)*bvec1(2,v) * b2*bvec2(1,v)^2), ... %xxxy
                   2 * (b1*bvec1(1,v)^2 * b2*bvec2(1,v)*bvec2(3,v) + ...
                        b1*bvec1(1,v)*bvec1(3,v) * b2*bvec2(1,v)^2), ... %xxxz
                   2 * (b1*bvec1(2,v)^2 * b2*bvec2(2,v)*bvec2(1,v) + ...
                        b1*bvec1(2,v)*bvec1(1,v) * b2*bvec2(2,v)^2), ... %xyyy
                   2 * (b1*bvec1(2,v)^2 * b2*bvec2(2,v)*bvec2(3,v) + ...
                        b1*bvec1(2,v)*bvec1(3,v) * b2*bvec2(2,v)^2), ... %yyyz (7)
                   2 * (b1*bvec1(3,v)^2 * b2*bvec2(3,v)*bvec2(1,v) + ...
                        b1*bvec1(3,v)*bvec1(1,v) * b2*bvec2(3,v)^2), ... %xzzz (8)
                   2 * (b1*bvec1(3,v)^2 * b2*bvec2(3,v)*bvec2(2,v) + ...
                        b1*bvec1(3,v)*bvec1(2,v) * b2*bvec2(3,v)^2), ... %zzzy (9)
                   b1*bvec1(1,v)^2 * b2*bvec2(2,v)^2 + ...
                       b1*bvec1(2,v)^2 * b2*bvec2(1,v)^2 , ... %xxyy (10)
                   b1*bvec1(1,v)^2 * b2*bvec2(3,v)^2 + ...
                       b1*bvec1(3,v)^2 * b2*bvec2(1,v)^2 , ... %xxzz (11)
                   b1*bvec1(2,v)^2 * b2*bvec2(3,v)^2 + ...
                       b1*bvec1(3,v)^2 * b2*bvec2(2,v)^2 , ... %yyzz (12)
                   2 * (b1*bvec1(1,v)^2 * b2*bvec2(2,v)*bvec2(3,v) + ...
                        b1*bvec1(2,v)*bvec1(3,v) * b2*bvec2(1,v)^2), ... %xxyz (13)
                   2 * (b1*bvec1(2,v)^2 * b2*bvec2(1,v)*bvec2(3,v) + ...
                        b1*bvec1(1,v)*bvec1(3,v) * b2*bvec2(2,v)^2), ... %yyxz (14)
                   2 * (b1*bvec1(3,v)^2 * b2*bvec2(1,v)*bvec2(2,v) + ...
                        b1*bvec1(1,v)*bvec1(2,v) * b2*bvec2(3,v)^2), ... %zzxy (15)
                   4 * (b1*bvec1(1,v)*bvec1(2,v) * b2*bvec2(1,v)*bvec2(2,v)),... %xyxy
                   4 * (b1*bvec1(1,v)*bvec1(3,v) * b2*bvec2(1,v)*bvec2(3,v)),... %xzxz
                   4 * (b1*bvec1(2,v)*bvec1(3,v) * b2*bvec2(2,v)*bvec2(3,v)),... %yzyz
                   4 * (b1*bvec1(1,v)*bvec1(2,v) * b2*bvec2(2,v)*bvec2(3,v) + ...
                        b1*bvec1(2,v)*bvec1(3,v) * b2*bvec2(1,v)*bvec2(2,v)) ... %xyyz
                   4 * (b1*bvec1(1,v)*bvec1(3,v) * b2*bvec2(2,v)*bvec2(3,v) + ...
                        b1*bvec1(2,v)*bvec1(3,v) * b2*bvec2(1,v)*bvec2(3,v)) ... %xzzy
                   4 * (b1*bvec1(1,v)*bvec1(2,v) * b2*bvec2(1,v)*bvec2(3,v) + ...
                        b1*bvec1(1,v)*bvec1(3,v) * b2*bvec2(1,v)*bvec2(2,v)) ... %yxxz
];

end

A=[-Ad Ak/6 Az ones(Nvol, 1)];





