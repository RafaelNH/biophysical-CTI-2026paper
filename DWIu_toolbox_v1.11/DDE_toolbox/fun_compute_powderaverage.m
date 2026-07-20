function [Data_pa] = fun_compute_powderaverage(data, gtab)

sized = size(data);
nb = length([0, gtab.ub1]);
sized(4) = nb;

bvals1_all = gtab.bval1;
bvals2_all = gtab.bval2;
dir1_all = gtab.dir1;
dir2_all = gtab.dir2;
b1_all = [0, gtab.ub1];
b2_all = [0, gtab.ub2];
cos2ang = [1, gtab.uc2a];

Data_pa_nonorm = zeros(sized);
for bi = 0:(nb-1)
    Data_pa_nonorm(:, :, :, bi+1) = mean(data(:, :, :, gtab.shells==bi), 4);
end

mS0 = mean(data(:, :, :, gtab.bval==0), 4);
for vi = 1:size(data, 4)
    data(:, :, :, vi) = data(:, :, :, vi) ./ mS0;
end
Data_all = data;


Data_pa = zeros(sized);
for bi = 0:(nb-1)
    Data_pa(:, :, :, bi+1) = mean(data(:, :, :, gtab.shells==bi), 4);
end

figure
for bi = 1:nb
    subplot(nb, 1, bi)
    data_si = Data_pa(:, :, :, bi);
    imagesc(data_si(:, :), [0 1.2])
end

dd = diag(dir1_all'*dir2_all);

figure, plot(bvals1_all), hold on, plot(bvals2_all, '--'), plot(dd - 2)