clear all
close all
clc

Ts = 2e-3;
z = tf('z', Ts);

G = zpk(minreal(0.02*(z+0.85)/((z-1) * (z-0.6))))
G1 = zpk(minreal((z-1) * G))
dcgain(G1)