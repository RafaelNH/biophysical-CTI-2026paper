pth = 'C:\Users\rafae\Data\Templates\CTIinvariants\';

v = load_untouch_nii([pth, 'MASK.nii']);
mask = double(v.img);


v = load_untouch_nii([pth, 'MD.nii']);
MD = double(v.img);
MD(mask==0) = 0;
v = load_untouch_nii([pth, 'AD.nii']);
AD = double(v.img);
AD(mask==0) = 0;
v = load_untouch_nii([pth, 'RD.nii']);
RD = double(v.img);
RD(mask==0) = 0;
v = load_untouch_nii([pth, 'FA.nii']);
FA = double(v.img);
FA(mask==0) = 0;

v = load_untouch_nii([pth, 'MK.nii']);
MK = double(v.img);
MK(mask==0) = 0;
v = load_untouch_nii([pth, 'AK.nii']);
AK = double(v.img);
AK(mask==0) = 0;
v = load_untouch_nii([pth, 'RK.nii']);
RK = double(v.img);
RK(mask==0) = 0;

v = load_untouch_nii([pth, 'uMK.nii']);
MKi = double(v.img);
MKi(mask==0) = 0;
v = load_untouch_nii([pth, 'uAK.nii']);
AKi = double(v.img);
AKi(mask==0) = 0;
v = load_untouch_nii([pth, 'uRK.nii']);
RKi = double(v.img);
RKi(mask==0) = 0;

v = load_untouch_nii([pth, 'vMK.nii']);
MKv = double(v.img);
MKv(mask==0) = 0;
v = load_untouch_nii([pth, 'vAK.nii']);
AKv = double(v.img);
AKv(mask==0) = 0;
v = load_untouch_nii([pth, 'vRK.nii']);
RKv = double(v.img);
RKv(mask==0) = 0;

v = load_untouch_nii([pth, 'KANI.nii']);
KANI = double(v.img);
KANI(mask==0) = 0;

v = load_untouch_nii([pth, 'KISO.nii']);
KISO = double(v.img);
KISO(mask==0) = 0;

v = load_untouch_nii([pth, 'OP.nii']);
OP = double(v.img);

v = load_untouch_nii([pth, 'uFA.nii']);
uFA = double(v.img);

v = load_untouch_nii('C:\Users\rafae\Data\Templates\JHU-ICBM-labels-2mm.nii');
ROI_temp = double(v.img);
