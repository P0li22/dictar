clear all
close all
clc

Ts = 0.05;
A = [-1.5 0; 1 0];
B = [1; 0];
C = [0 1];
D = 0;
x0 = [0; 0];
sys = ss(A, B, C, D);
sys_x = ss(A, B, eye(2), 0);
sys_d = c2d(sys, Ts, 'zoh');
Ad = sys_d.a;
Bd = sys_d.b;
Cd = sys_d.c;
Dd = sys_d.d;

Cy = eye(2);
Cz = Cd;
Cc = Cz;
Dz = Dd;
Dc = Dz;

Q = 10;
R = 1;

Hw = 1;
Hp = 20;
Hc = 20;

zblk = 1;
ublk = 1;
W = []; 
V = [];
cmode=0;
solver = 'qp_as';

du_max = [inf]; du_min = [-inf];
u_max = [1]; u_min = [-1];
z_max = [1.0099]; z_min = [0];

md = MPCInit(Ad, Bd, Cy, Cz, Dz, Cc, Dc, Hp, Hw, zblk, Hc, ublk, du_max, du_min, ...
    u_max, u_min, z_max, z_min, Q, R, W, V, Ts, cmode, solver);

out=sim('MPC_model');
figure(1)
subplot(221)
plot(y.time, y.data);
hold on, grid on, zoom on,
xline(2.5, ':r');
yline(1.05, ':r');
yline(1.01, ':r'); yline(0.99, ':r');
xlabel('t'); ylabel('y');

subplot(222)
plot(u.time, u.data);
hold on, grid on, zoom on, 
yline(1, ':r');
xlabel('t'); ylabel('u');



















