function bvec1 = fun_bvec1_from_bvec2_datatype1(bvec2)

% function to reconstruct bvec1 from bvec2 given that data contains the
% following sequence of volumes:
%     a) 1 b=0 s/mm^2
%     b) 61 b1=2000 s/mm^2 b2=0    s/mm^2 corresponding to set #2 in Fig. 1 NeuroImage paper (including b=0 volumes at the beginning or at the end), POSITIVE POLARITY
%     c) 61 b1=1000 s/mm^2 b2=1000 s/mm^2 corresponding to set #3 in Fig. 1 NeuroImage paper (including b=0 volumes at the beginning or at the end) parallel directions, POSITIVE POLARITY
%     d) 61 b1=1000 s/mm^2 b2=1000 s/mm^2 corresponding to set #4 in Fig. 1 NeuroImage paper (including b=0 volumes at the beginning or at the end) perpendicular directions, POSITIVE POLARITY
%     e) 61 b1=1000 s/mm^2 b2=0    s/mm^2 corresponding to set #1 in Fig. 1 NeuroImage paper (including b=0 volumes at the beginning or at the end), POSITIVE POLARITY
%     f) 61 b1=2000 s/mm^2 b2=0    s/mm^2 corresponding to set #2 in Fig. 1 NeuroImage paper (including b=0 volumes at the beginning or at the end), NEGATIVE POLARITY
%     g) 61 b1=1000 s/mm^2 b2=1000 s/mm^2 corresponding to set #3 in Fig. 1 NeuroImage paper (including b=0 volumes at the beginning or at the end) parallel directions, NEGATIVE POLARITY
%     h) 61 b1=1000 s/mm^2 b2=1000 s/mm^2 corresponding to set #4 in Fig. 1 NeuroImage paper (including b=0 volumes at the beginning or at the end) perpendicular directions, NEGATIVE POLARITY
%     i) 61 b1=1000 s/mm^2 b2=0    s/mm^2 corresponding to set #1 in Fig. 1 NeuroImage paper (including b=0 volumes at the beginning or at the end), NEGATIVE POLARITY

% for parallel experiments bvecs1 = bvecs2
bvec1 = bvec2;

% perpendicular experiments (2 sets)
vvv_set1 = bvec2(:, 125:184)';
vvv_set2 = bvec2(:, 369:428)';

% reference directions
[V1, V2] = DDE_5design;
V1(:, 3) = -V1(:, 3);
V2(:, 3) = -V2(:, 3);

% rotate q1 directiions according to the rotation matrix between q2
% directions and its reference
bvecs1_sel = zeros(size(vvv_set1));
for vi=1:60
    R = findRotationMatrix(V2(vi, :), vvv_set1(vi, :));
    v1_ref_sel = V1(vi, :);
    v1_ref_sel = v1_ref_sel / norm(v1_ref_sel);
    bvecs1_sel(vi, :) = v1_ref_sel*R;
end

bvec1(:, 125:184) = bvecs1_sel';

for vi=1:60
    R = findRotationMatrix(V2(vi, :), vvv_set2(vi, :));
    v1_ref_sel = V1(vi, :);
    v1_ref_sel = v1_ref_sel / norm(v1_ref_sel);
    bvecs1_sel(vi, :) = v1_ref_sel*R;
end

bvec1(:, 369:428) = bvecs1_sel';

