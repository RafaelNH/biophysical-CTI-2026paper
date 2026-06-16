function [grad1, grad2] = dde_profile(g1, g2, delta, Delta, rsamp)

len2 = (delta+Delta*3+2*delta)*rsamp;
len1 = (delta+Delta*3+2*delta)*rsamp/2;

grad2 = zeros(1, len2);
grad1 = zeros(1, len1);


grad1(1, (delta*rsamp+1):(2*delta)*rsamp) = g1;
grad1(1, (delta*rsamp+Delta*rsamp+1):(2*delta*rsamp+Delta*rsamp)) = -g1;

grad2(1, (2*Delta*rsamp+delta*rsamp+1):(2*Delta*rsamp+2*delta*rsamp)) = g2;
grad2(1, (delta*rsamp+3*Delta*rsamp+1):(2*delta*rsamp+3*Delta*rsamp)) = -g2;
grad2(1, 1:len1-1) = nan;
end