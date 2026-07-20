function gtab = fun_reconst_gtab_for_human(btotal, bvecs2, datatype)
%%
%gtab.ub1 %b1s of each shell
%gtab.ub2 %b2s of each shell
%gtab.uc2a % cos2ang of each shell
%
%gtab.bval1 % all b1s
%gtab.bval2 % all b2s
%gtab.bval % all total bvalues
%gtab.dir1 % all directions 1
%gtab.dir2 % all directions 2
%gtab.B1 % bmatrix 1 of individual data
%gtab.B2 % bmatrix 2 of individual data
%gtab.c2a % cos2ang of each image
%
%gtab.shells % varianble indicating the indexes
% 0 is b1=b2=0, indexs 1 to 4 corresponds to the 4 sets of DDE experiments
%%

if datatype == 1
    [bvals1, bvals2, ub1, ub2, uc2a, shells] = fun_bt_2_b1b2_datatype1(btotal);
elseif datatype == 2
    [bvals1, bvals2, ub1, ub2, uc2a, shells] = fun_bt_2_b1b2_datatype2(btotal);
end

bvecs1 = fun_bvec1_from_bvec2_v2(bvecs2, shells);

gtab.uc2a = uc2a;
gtab.ub1 = ub1;
gtab.ub2 = ub2;

gtab.bval1 = bvals1; % all b1s
gtab.bval2 = bvals2; % all b2s
gtab.bval = btotal; % all total bvalues
gtab.dir1 = bvecs1;  % all directions 1
gtab.dir2 = bvecs2;  % all directions 2
gtab.c2a = diag(bvecs1' * bvecs2)';
gtab.shells = shells;