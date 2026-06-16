function [dt, kt] = fun_rotate_tensors_DKI(dt, kt, Evector)

DTensor = [dt(1), dt(4), dt(5); ...
    dt(4), dt(2), dt(6);...
    dt(5), dt(6), dt(3)];
DTr = Evector * DTensor * pinv(Evector);
dt = [DTr(1, 1), DTr(2, 2), DTr(3, 3), DTr(1, 2), DTr(1, 3), DTr(2, 3)];

W1111_r = Wrotate(kt, 1, 1, 1, 1, Evector');
W2222_r = Wrotate(kt, 2, 2, 2, 2, Evector');
W3333_r = Wrotate(kt, 3, 3, 3, 3, Evector');
W1112_r = Wrotate(kt, 1, 1, 1, 2, Evector');
W1113_r = Wrotate(kt, 1, 1, 1, 3, Evector');
W1222_r = Wrotate(kt, 1, 2, 2, 2, Evector');
W2223_r = Wrotate(kt, 2, 2, 2, 3, Evector');
W1333_r = Wrotate(kt, 1, 3, 3, 3, Evector');
W2333_r = Wrotate(kt, 2, 3, 3, 3, Evector');
W1122_r = Wrotate(kt, 1, 1, 2, 2, Evector');
W1133_r = Wrotate(kt, 1, 1, 3, 3, Evector');
W2233_r = Wrotate(kt, 2, 2, 3, 3, Evector');
W1123_r = Wrotate(kt, 1, 1, 2, 3, Evector');
W1223_r = Wrotate(kt, 1, 2, 2, 3, Evector');
W1233_r = Wrotate(kt, 1, 2, 3, 3, Evector');

kt = [W1111_r, W2222_r, W3333_r, W1112_r, W1113_r, W1222_r,...
    W2223_r, W1333_r, W2333_r, W1122_r, W1133_r, W2233_r, W1123_r, W1223_r, W1233_r];

