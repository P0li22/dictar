clear all 
close all
clc

s = tf('s');
Gcont = 40 / ( s^2 + 4*s - 9.81 );

Ts = 0.2;
z = tf('z', Ts);
Gz = zpk( c2d(Gcont, Ts, 'zoh') )
zeta = 0.456;
%zeta = 0.6;
wn = 2.24;

% pzmap(Gz);
% hold on;
% zgrid(zeta, 0, Ts);

p1c = -zeta * wn + j * wn * sqrt(1 - zeta^2)
p2c = -zeta * wn - j * wn * sqrt(1 - zeta^2)
p3c = -5 * zeta * wn
p4c = -5 * zeta * wn

p1_d = exp(p1c * Ts)
p2_d = exp(p2c * Ts)
p3_d = exp(p3c * Ts)
p4_d = exp(p4c * Ts)

[zg, pg, kg] = zpkdata(Gz, 'v')
A = conv([1 -pg(1)], [1 -pg(2)]);
B = kg * [1 -zg];

A_plus = [1 -pg(2)]; 
A_minus = [1 -pg(1)];
B_plus = 1;
B_minus = B;
A_dioph = conv([1, -1], A_minus);
B_dioph = B_minus;
Am = poly([p1_d, p2_d, p3_d, p4_d]);
Ms = [ [A_dioph(:); 0; 0], [0; A_dioph(:); 0], [0; 0; A_dioph(:)]...
    [0; B_dioph(:); 0; 0], [0; 0; B_dioph(:); 0], [0; 0; 0; B_dioph(:)];...
    [0.03, 0.03, 0.03, -0.006812, -0.006812, -0.006812] ]
gamma = [Am'; 0]
theta = Ms\gamma

R1 = theta(1:3)'
S1 = theta(4:6)'
S = conv(A_plus, S1)
R = conv([1, -1], R1)
C = zpk(tf(S, R, Ts))



Wprime = (C*Gz)/(1+(C*Gz));

T_tilde = conv([1 -p3_d], [1 -p4_d]);
kW = 1;
kT = kW*polyval(S1,1)/(polyval(T_tilde,1)*dcgain(Wprime));

T1 = kT*T_tilde;
F = zpk(tf(T1,S1,Ts))

% performance evaluation
out = sim("x2dof_sim.slx");
figure(1)
title("y")
hold on
plot(out.y.Time, out.y.Data)
grid on
figure(2)
plot(out.u.Time, out.u.Data)
grid on























