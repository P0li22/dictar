clear all 
close all
clc

A = [0 1; 0 -1];
B = [0; 1];
C = [1 0];
D = 0;
sys = ss(A, B, C, D);
x0 = [0; 0];
Ts = 0.05;
sys_d = c2d(sys, Ts, 'zoh');

Ad = sys_d.A;
Bd = sys_d.B;
Cd = sys_d.C;
Dd = sys_d.D;
n = 2;

Mr = ctrb(Ad, Bd);
rank_mr = rank(Mr)

Aaug = [1 -Ts*Cd; zeros(n, 1) Ad];
Baug = [0; Bd];
Caug = [0 Cd];
Daug = 0;

Q = diag([10600 250 1]);
R = 1;

Cq = chol(Q);
Mo = obsv(Aaug, Cq);
rank_mo = rank(Mo)

Kaug = dlqr(Aaug, Baug, Q, R)
Kq = Kaug(1)
Kx = Kaug(2:n+1)

sys_x = ss(A, B, eye(n), 0);
out = sim('lab5_es1_sim')
figure(1),
plot(out.r.time, out.r.data)
grid on, zoom on, hold on,
plot(out.y.time, out.y.data)
xlabel('t'), ylabel('u and y')
xline(2, ':r')

figure(2)
subplot(221)
plot(out.x.time, out.x.data(:, 1))
grid on, zoom on, hold on,
xlabel('t'), ylabel('x1'),
subplot(222)
plot(out.x.time, out.x.data(:, 2))
grid on, zoom on, hold on, 
xlabel('t'), ylabel('x2'),
subplot(223)
plot(out.x.time, sqrt( out.x.data(:, 1).^2+out.x.data(:,2).^2 ))
grid on, zoom on, hold on
xlabel('t'), ylabel('||x||')
subplot(224)
plot(out.u.time, out.u.data)
grid on, zoom on, hold on,
xlabel('t'), ylabel('u')
yline(1, ':r')



















