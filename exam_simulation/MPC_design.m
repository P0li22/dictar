clear all
close all
clc

A = [0 1; 0 -0.1];
B = [0; 0.1];
C = [1 0];
D = 0;
x0 = [0; 0];
Ts = 50e-3;

sys = ss(A, B, C, D);
sys_x = ss(A, B, eye(2), 0);
sys_d = c2d(sys, Ts, 'zoh');

Ad = sys_d.a;
Bd = sys_d.b;
Cd = sys_d.c;
Dd = sys_d.d;

Cy = eye(2);
Cz = Cd;
Dz = 0;
Cc = Cz;
Dc = Dz;

Q = 500;
R = 1;

Hp = 35;
Hc = 10;
Hw = 1;
ublk = 1;
zblk = 1;
du_max = [2];
du_min = [-2];
u_max = [10];
u_min = [-10];
z_max = [1.0199];
z_min = [-inf];

md = MPCInit(Ad, Bd, Cy, Cz, Dz, Cc, Dc, Hp, Hw, zblk, Hc, ublk, du_max, ...
    du_min, u_max, u_min, z_max, z_min, Q, R, [], [], Ts, 0, 'qp_as')

out = sim("MPC_design_model.slx")
figure(1)
plot(out.y.time, out.y.data)
grid on, zoom on, hold on
xlabel('time'), ylabel('y')
yline(1+1*2/100, ':r'), yline(1-1*2/100, ':r'),
xline(2, ':r')

figure(2)
plot(out.u.time, out.u.data)
grid on, zoom on, hold on
xlabel('time'), ylabel('u')


