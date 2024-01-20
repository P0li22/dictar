clear all
close all
clc

A = [0 1; 0 -1];
B = [0; 1];
C = [1 0];
D = 0;
sys = ss(A, B, C, D);
x0 = [0.8; 0];
Ts = 0.05;
n = 2;
sys_dt = c2d(sys, Ts, 'zoh');

Ad = sys_dt.a;
Bd = sys_dt.b;
Cd = sys_dt.c;
Dd = sys_dt.d;

Mr = ctrb(Ad, Bd);
rho_Mr = rank(Mr)

Q = [300 0; 0 1];
R = 1;

Cq = chol(Q);
Mo = obsv(Ad, Cq);
rho_Mo = rank(Mo)
K = dlqr(Ad, Bd, Q, R);



sys_x = ss(A,B,eye(2),0);
t_sim = 3;
out = sim('lab4_sim');
figure(1);
subplot(221);
plot(out.x.time,out.x.data(:,1),'b','linew',1.2);
grid on, zoom on, hold on,
xlabel('t'),ylabel('x_1(t)')
subplot(222)
plot(out.x.time,out.x.data(:,2),'b','linew',1.2)
grid on,zoom on, hold on,
xlabel('t'),ylabel('x_2(t)')
subplot(223),
plot(out.x.time,sqrt(out.x.data(:,1).^2+out.x.data(:,2).^2),'b','linew',1.2)
hold on,zoom on, grid on
yline(1e-5,':r','linew',1.2)
yline(-1e-5,':r','linew',1.2)
xlabel('t'),ylabel('||x(t)||_2')
subplot(224),
plot(out.u.time,out.u.data,'b','linew',1.2)
grid on ,zoom on, hold on,
xlabel('t'),ylabel('u(t)')





