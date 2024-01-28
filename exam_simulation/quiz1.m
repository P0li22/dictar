clear all
close all 
clc

Ts = 1;
z = tf('z', Ts)
C = 1.5 * (z-0.7) / (z-0.9);
G = 1.1 * 1/(z-1.2);

L = C * G;
W = zpk(minreal(L / (1+L)));
W1 = zpk(minreal(G / (1+L)));
W2 = zpk(minreal(1 / (1+L)));

[zg, pg, kg] = zpkdata(W, 'v');

R = z / (z-1);
D1 = R;

Y = zpk(minreal(W * R + W1 * D1))
[num, den] = tfdata(Y, 'v');
[r, p, k] = residuez(num, den)