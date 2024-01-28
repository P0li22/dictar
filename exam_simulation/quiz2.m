clear all
close all
clc

Ts = 1;
z = tf('z', Ts);
C = 0.4 * (z-0.4) / (z-0.2);
G = (z+0.2) / ( (z-1) * (z-0.1) );
L = C * G;
W = zpk(minreal(L/(1+L)))
W1 = zpk(minreal(G/(1+L)))
W2 = zpk(minreal(1/(1+L)))

[zg, pg, kg] = zpkdata(W, 'v')
R = 0;
D1 = 0.15 * z/(z-1);
D2 = 0.1 * z / (z-1);

Y = zpk(minreal(W * R + W1 * D1 + W2 * D2));
Y1 = zpk(minreal((z-1) * Y));

yinf = dcgain(Y1)