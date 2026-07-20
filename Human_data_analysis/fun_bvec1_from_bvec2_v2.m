function bvec1 = fun_bvec1_from_bvec2_v2(bvec2, shells)

% function to reconstruct bvec1 from bvec2 assuming that shells 3 and 7
% corresponds to perpendicular DDE experiments.

% Normalize bvecs
for vi = 1:size(bvec2, 2)
    if norm(bvec2(:, vi))> 0
        bvec2(:, vi) = bvec2(:, vi) / norm(bvec2(:, vi));
    end
end

bvec1 = bvec2;

% perpendicular experiments (2 sets)
vvv_set1 = bvec2(:, shells==3)';
vvv_set2 = bvec2(:, shells==7)';

% reference directions
[V1, V2] = DDE_5design;
V1(:, 3) = -V1(:, 3); % scanner applies negative direction in z
V2(:, 3) = -V2(:, 3); % scanner applies negative direction in z
V2(37, 3) = -V2(37, 3); % There was a typo and this gradient was acquired
% with inverted z element

figure('color', [1 1 1]),
plot(acos(diag(V2*vvv_set1'))/pi*180, 'blue')
hold on
plot(acos(diag((-V2)*vvv_set2'))/pi*180, 'red')
yline(1)

% % rotate q1 directiions according to the rotation matrix between q2
% % directions and its reference
% bvecs1_sel = zeros(size(vvv_set1));
% for vi=1:60
%     R = findRotationMatrix(V2(vi, :), vvv_set1(vi, :));
%     v1_ref_sel = V1(vi, :);
%     v1_ref_sel = v1_ref_sel / norm(v1_ref_sel);
%     bvecs1_sel(vi, :) = v1_ref_sel*R;
% end

%bvec1(:, shells==3) = bvecs1_sel';
bvec1(:, shells==3) = V1';

% for vi=1:60
%     R = findRotationMatrix(V2(vi, :), vvv_set2(vi, :));
%     v1_ref_sel = V1(vi, :);
%     v1_ref_sel = v1_ref_sel / norm(v1_ref_sel);
%     bvecs1_sel(vi, :) = v1_ref_sel*R;
% end

% figure, plot(diag(V1*bvecs1_sel'))

%bvec1(:, shells==7) = bvecs1_sel';
bvec1(:, shells==7) = -V1';

