clear all
close all
clc

z = tf('z', 1);
C = 0.2*(z+0.1)/(z-0.2);
G = 1/(z-1);

D1 = 2*z/(z-1);
L = zpk(minreal(C*G));
W1 = zpk(minreal(G/(1+L)));
Y = zpk(minreal(W1 * D1));
[num, den] = tfdata(Y, 'v');
[r, p, k] = residuez(num, den)