function dij = fun_reconst_d_from_pos(R, time, i, j)

aRiRj = mean(R(:, i) .* R(:, j));
dij = aRiRj / (2*time);