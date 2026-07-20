function [b1, b2, ub1, ub2, uc2a, shells] = fun_bt_2_b1b2_datatype2(btotal)
% function to convert total DDE bval to b1 / b2
% given that data contains the following sequence of volumes:
%a) 1 b=0 s/mm^2
%b) 66 b1=2000 s/mm^2 b2=0    s/mm^2 corresponding to set #2 in Fig. 1 NeuroImage paper (including b=0 volumes interspersed), POSITIVE POLARITY
%c) 65 b1=1000 s/mm^2 b2=1000 s/mm^2 corresponding to set #3 in Fig. 1 NeuroImage paper (including b=0 volumes interspersed) parallel directions, POSITIVE POLARITY
%d) 66 b1=1000 s/mm^2 b2=1000 s/mm^2 corresponding to set #4 in Fig. 1 NeuroImage paper (including b=0 volumes interspersed) perpendicular directions, POSITIVE POLARITY
%e) 65 b1=1000 s/mm^2 b2=0    s/mm^2 corresponding to set #1 in Fig. 1 NeuroImage paper (including b=0 volumes interspersed), POSITIVE POLARITY
%f) 66 b1=2000 s/mm^2 b2=0    s/mm^2 corresponding to set #2 in Fig. 1 NeuroImage paper (including b=0 volumes interspersed), NEGATIVE POLARITY
%g) 65 b1=1000 s/mm^2 b2=1000 s/mm^2 corresponding to set #3 in Fig. 1 NeuroImage paper (including b=0 volumes interspersed) parallel directions, NEGATIVE POLARITY
%h) 66 b1=1000 s/mm^2 b2=1000 s/mm^2 corresponding to set #4 in Fig. 1 NeuroImage paper (including b=0 volumes interspersed) perpendicular directions, NEGATIVE POLARITY
%i) 65 b1=1000 s/mm^2 b2=0    s/mm^2 corresponding to set #1 in Fig. 1 NeuroImage paper (including b=0 volumes interspersed), NEGATIVE POLARITY
% 

% variable wset defines which type of experiment correspond each volume
% based on the information above
wset = zeros(size(btotal));
shells = zeros(size(btotal));

% wset = 0 for b0s
% wset = 1 for b1>0 and b2=0 experiments
% wset = 2 for b1 = b2 experiments
wset(3:67) = 1; %b)
wset(68:132) = 2; %c)
wset(133:198) = 2; %d) perpendicular
wset(199:263) = 1; %e)

wset(264:329) = 1; %f)
wset(330:394) = 2; %g)
wset(395:460) = 2; %h)
wset(461:525) = 1; %i)

shells(3:66) = 1; %b)
shells(68:132) = 2; %c)
shells(133:198) = 3; %d) perpendicular
shells(199:263) = 4; %e)

shells(264:329) = 5; %f)
shells(330:394) = 6; %g)
shells(395:460) = 7; %h) perpendicular
shells(461:525) = 8; %i)

% now were are ready to define b1 and b2 according to the different cases
b1 = zeros(size(btotal));
b2 = zeros(size(btotal));

b1(wset==1) = btotal(wset==1);
b1(wset==2) = btotal(wset==2)/2;
b2(wset==2) = btotal(wset==2)/2;

shells(btotal==0) = 0;

ub1 = b1([3, 68, 134, 199, 265, 330, 396, 461]);
ub2 = b2([3, 68, 134, 199, 265, 330, 396, 461]);
uc2a = [1, 1, 0, 1, 1, 1, 0, 1];

figure('color', [1 1 1])
plot(btotal, 'g')

hold on
plot(b1, 'b')
plot(b2, '.r')
plot(shells * 100, 'black')

