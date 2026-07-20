function fun_save_nii(Metric, vnii, save_name)

vnii.hdr.dime.dim(5)=1;
vnii.hdr.dime.pixdim(5)=0;
vnii.hdr.dime.xyz_t=0;
vnii.hdr.dime.datatype=16;
vnii.img=Metric;
save_untouch_nii(vnii, save_name)