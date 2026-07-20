function [b1, b2, ub1, ub2, uc2a, shells] = fun_bt_2_b1b2_datatype1(btotal)
% function to convert total DDE bval to b1 / b2
% given that data contains the following sequence of volumes:
%     a) 1 b=0 s/mm^2
%     b) 61 b1=2000 s/mm^2 b2=0    s/mm^2 corresponding to set #2 in Fig. 1 NeuroImage paper (including b=0 volumes at the beginning or at the end), POSITIVE POLARITY
%     c) 61 b1=1000 s/mm^2 b2=1000 s/mm^2 corresponding to set #3 in Fig. 1 NeuroImage paper (including b=0 volumes at the beginning or at the end) parallel directions, POSITIVE POLARITY
%     d) 61 b1=1000 s/mm^2 b2=1000 s/mm^2 corresponding to set #4 in Fig. 1 NeuroImage paper (including b=0 volumes at the beginning or at the end) perpendicular directions, POSITIVE POLARITY
%     e) 61 b1=1000 s/mm^2 b2=0    s/mm^2 corresponding to set #1 in Fig. 1 NeuroImage paper (including b=0 volumes at the beginning or at the end), POSITIVE POLARITY
%     f) 61 b1=2000 s/mm^2 b2=0    s/mm^2 corresponding to set #2 in Fig. 1 NeuroImage paper (including b=0 volumes at the beginning or at the end), NEGATIVE POLARITY
%     g) 61 b1=1000 s/mm^2 b2=1000 s/mm^2 corresponding to set #3 in Fig. 1 NeuroImage paper (including b=0 volumes at the beginning or at the end) parallel directions, NEGATIVE POLARITY
%     h) 61 b1=1000 s/mm^2 b2=1000 s/mm^2 corresponding to set #4 in Fig. 1 NeuroImage paper (including b=0 volumes at the beginning or at the end) perpendicular directions, NEGATIVE POLARITY
%     i) 61 b1=1000 s/mm^2 b2=0    s/mm^2 corresponding to set #1 in Fig. 1 NeuroImage paper (including b=0 volumes at the beginning or at the end), NEGATIVE POLARITY
%
% 

% variable wset defines which type of experiment correspond each volume
% based on the information above
wset = zeros(size(btotal));
shells = zeros(size(btotal));

% wset = 0 for b0s
% wset = 1 for b1>0 and b2=0 experiments
% wset = 2 for b1 = b2 experiments
wset(3:62) = 1; %b)
wset(63:122) = 2; %c)
wset(125:184) = 2; %d) perpendicular
wset(185:244) = 1; %d)

wset(247:306) = 1; %b)
wset(307:366) = 2; %c)
wset(369:428) = 2; %d)
wset(429:488) = 1; %d)

shells(3:62) = 1; %b)
shells(63:122) = 2; %c)
shells(125:184) = 3; %d) perpendicular
shells(185:244) = 4; %d)

shells(247:306) = 5; %b)
shells(307:366) = 6; %c)
shells(369:428) = 7; %d)
shells(429:488) = 8; %d)

% now were are ready to define b1 and b2 according to the different cases
b1 = zeros(size(btotal));
b2 = zeros(size(btotal));

b1(wset==1) = btotal(wset==1);
b1(wset==2) = btotal(wset==2)/2;
b2(wset==2) = btotal(wset==2)/2;

ub1 = b1([3, 63, 125, 185, 247, 307, 369, 429]);
ub2 = b2([3, 63, 125, 185, 247, 307, 369, 429]);
uc2a = [1, 1, 0, 1, 1, 1, 0, 1];

figure('color', [1 1 1])
plot(btotal, 'g')

hold on
plot(b1, 'b')
plot(b2, '.r')
plot(shells * 100, 'black')

