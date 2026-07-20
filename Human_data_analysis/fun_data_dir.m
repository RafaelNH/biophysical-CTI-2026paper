function [path, dwi_name, bval_name, bvec_name, mask_name, datatype] = fun_data_dir(whichd)

if whichd == 1
    path = 'C:\Users\rafae\Data\CTIhuman\sub-01\';
    dwi_name = 'sub-01\sub-01_dde-dwi_den_grc_tec-preproc-3_all-ser_all-sh_pos-neg-conc_drift-corr_bias-corr.nii';
    mask_name = 'sub-01\sub-01_dde-dwi_den_grc_tec-preproc-3_all-ser_all-sh_pos-neg-conc_drift-corr_bias-corr_brain_mask.nii';
    bval_name = 'sub-01_dde-dwi_den_grc_tec-preproc-3_all-ser_all-sh_pos-neg-conc_mean-bvals.txt';
    bvec_name = 'sub-01_dde-dwi_den_grc_tec-preproc-3_all-ser_all-sh_pos-neg-conc_gt-ed-rot-bvecs.txt';
    datatype = 1;
elseif whichd == 2
    path = 'C:\Users\rafae\Data\CTIhuman\sub-02\';
    dwi_name = 'sub-02\sub-02_dde-dwi_den_grc_tec-preproc-3_all-ser_all-sh_pos-neg-conc_drift-corr_bias-corr.nii';
    mask_name = 'sub-02\sub-02_dde-dwi_den_grc_tec-preproc-3_all-ser_all-sh_pos-neg-conc_drift-corr_bias-corr_brain_mask.nii';
    bval_name = 'sub-02_dde-dwi_den_grc_tec-preproc-3_all-ser_all-sh_pos-neg-conc_mean-bvals.txt';
    bvec_name = 'sub-02_dde-dwi_den_grc_tec-preproc-3_all-ser_all-sh_pos-neg-conc_gt-ed-rot-bvecs.txt';
    datatype = 1;
elseif whichd == 3
    path = 'C:\Users\rafae\Data\CTIhuman\sub-03\';
    dwi_name = 'sub-03\sub-03_dde-dwi_den_grc_tec-preproc-3_all-ser_all-sh_pos-neg-conc_drift-corr_bias-corr.nii';
    mask_name = 'sub-03\sub-03_dde-dwi_den_grc_tec-preproc-3_all-ser_all-sh_pos-neg-conc_drift-corr_bias-corr_brain_mask.nii';
    bval_name = 'sub-03_dde-dwi_den_grc_tec-preproc-3_all-ser_all-sh_pos-neg-conc_mean-bvals.txt';
    bvec_name = 'sub-03_dde-dwi_den_grc_tec-preproc-3_all-ser_all-sh_pos-neg-conc_gt-ed-rot-bvecs.txt';
    datatype = 2;
elseif whichd == 4
    path = 'C:\Users\rafae\Data\CTIhuman\sub-04\';
    dwi_name = 'sub-04\sub-04_dde-dwi_den_grc_tec-preproc-3_all-ser_all-sh_pos-neg-conc_drift-corr_bias-corr.nii';
     mask_name = 'sub-04\sub-04_dde-dwi_den_grc_tec-preproc-3_all-ser_all-sh_pos-neg-conc_drift-corr_bias-corr_brain_mask.nii';
    bval_name = 'sub-04_dde-dwi_den_grc_tec-preproc-3_all-ser_all-sh_pos-neg-conc_mean-bvals.txt';
    bvec_name = 'sub-04_dde-dwi_den_grc_tec-preproc-3_all-ser_all-sh_pos-neg-conc_gt-ed-rot-bvecs.txt';
    datatype = 2;
elseif whichd == 5
    path = 'C:\Users\rafae\Data\CTIhuman\sub-05\';
    dwi_name = 'sub-05\sub-05_dde-dwi_den_grc_tec-preproc-3_all-ser_all-sh_pos-neg-conc_drift-corr_bias-corr.nii';
    mask_name = 'sub-05\sub-05_dde-dwi_den_grc_tec-preproc-3_all-ser_all-sh_pos-neg-conc_drift-corr_bias-corr_brain_mask.nii';
    bval_name = 'sub-05_dde-dwi_den_grc_tec-preproc-3_all-ser_all-sh_pos-neg-conc_mean-bvals.txt';
    bvec_name = 'sub-05_dde-dwi_den_grc_tec-preproc-3_all-ser_all-sh_pos-neg-conc_gt-ed-rot-bvecs.txt';
    datatype = 2;
elseif whichd == 6
    path = 'C:\Users\rafae\Data\CTIhuman\sub-06\';
    dwi_name = 'sub-06\sub-06_dde-dwi_den_grc_tec-preproc-3_all-ser_all-sh_pos-neg-conc_drift-corr_bias-corr.nii';
    mask_name = 'sub-06\sub-06_dde-dwi_den_grc_tec-preproc-3_all-ser_all-sh_pos-neg-conc_drift-corr_bias-corr_brain_mask.nii';
    bval_name = 'sub-06_dde-dwi_den_grc_tec-preproc-3_all-ser_all-sh_pos-neg-conc_mean-bvals.txt';
    bvec_name = 'sub-06_dde-dwi_den_grc_tec-preproc-3_all-ser_all-sh_pos-neg-conc_gt-ed-rot-bvecs.txt';
    datatype = 2;
elseif whichd == 7
    path = 'C:\Users\rafae\Data\CTIhuman\sub-07\';
    dwi_name = 'sub-07\sub-07_dde-dwi_den_grc_tec-preproc-3_all-ser_all-sh_pos-neg-conc_drift-corr_bias-corr.nii';
    mask_name = 'sub-07\sub-07_dde-dwi_den_grc_tec-preproc-3_all-ser_all-sh_pos-neg-conc_drift-corr_bias-corr_brain_mask.nii';
    bval_name = 'sub-07_dde-dwi_den_grc_tec-preproc-3_all-ser_all-sh_pos-neg-conc_mean-bvals.txt';
    bvec_name = 'sub-07_dde-dwi_den_grc_tec-preproc-3_all-ser_all-sh_pos-neg-conc_gt-ed-rot-bvecs.txt';
    datatype = 2;
elseif whichd == 8
    path = 'C:\Users\rafae\Data\CTIhuman\sub-08\';
    dwi_name = 'sub-08\sub-08_dde-dwi_den_grc_tec-preproc-3_all-ser_all-sh_pos-neg-conc_drift-corr_bias-corr.nii';
    mask_name = 'sub-08\sub-08_dde-dwi_den_grc_tec-preproc-3_all-ser_all-sh_pos-neg-conc_drift-corr_bias-corr_brain_mask.nii';
    bval_name = 'sub-08_dde-dwi_den_grc_tec-preproc-3_all-ser_all-sh_pos-neg-conc_mean-bvals.txt';
    bvec_name = 'sub-08_dde-dwi_den_grc_tec-preproc-3_all-ser_all-sh_pos-neg-conc_gt-ed-rot-bvecs.txt';
    datatype = 2;
elseif whichd == 9
    path = 'C:\Users\rafae\Data\CTIhuman\sub-09\';
    dwi_name = 'sub-09\sub-09_dde-dwi_den_grc_tec-preproc-3_all-ser_all-sh_pos-neg-conc_drift-corr_bias-corr.nii';
    mask_name = 'sub-09\sub-09_dde-dwi_den_grc_tec-preproc-3_all-ser_all-sh_pos-neg-conc_drift-corr_bias-corr_brain_mask.nii';
    bval_name = 'sub-09_dde-dwi_den_grc_tec-preproc-3_all-ser_all-sh_pos-neg-conc_mean-bvals.txt';
    bvec_name = 'sub-09_dde-dwi_den_grc_tec-preproc-3_all-ser_all-sh_pos-neg-conc_gt-ed-rot-bvecs.txt';
    datatype = 2;
elseif whichd == 10
    path = 'C:\Users\rafae\Data\CTIhuman\sub-10\';
    dwi_name = 'sub-10\sub-10_dde-dwi_den_grc_tec-preproc-3_all-ser_all-sh_pos-neg-conc_drift-corr_bias-corr.nii';
    mask_name = 'sub-10\sub-10_dde-dwi_den_grc_tec-preproc-3_all-ser_all-sh_pos-neg-conc_drift-corr_bias-corr_brain_mask.nii';
    bval_name = 'sub-10_dde-dwi_den_grc_tec-preproc-3_all-ser_all-sh_pos-neg-conc_mean-bvals.txt';
    bvec_name = 'sub-10_dde-dwi_den_grc_tec-preproc-3_all-ser_all-sh_pos-neg-conc_gt-ed-rot-bvecs.txt';
    datatype = 2;
end

    