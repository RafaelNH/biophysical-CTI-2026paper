function kijkl = fun_reconst_k_from_pos(R, ii, jj, kk, ll)


aRiRjRkRl = mean(R(:, ii) .* R(:, jj) .* R(:, kk) .* R(:, ll));
aRiRj = mean(R(:, ii) .* R(:, jj));
aRkRl = mean(R(:, kk) .* R(:, ll));
aRiRk = mean(R(:, ii) .* R(:, kk));
aRjRl = mean(R(:, jj) .* R(:, ll));
aRiRl = mean(R(:, ii) .* R(:, ll));
aRjRk = mean(R(:, jj) .* R(:, kk));

kijkl = aRiRjRkRl - aRiRj * aRkRl ...
                  - aRiRk * aRjRl ...
                  - aRiRl * aRjRk;