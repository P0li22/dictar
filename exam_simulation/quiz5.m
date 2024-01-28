clear all
close all
clc

A = 0.732;
B = 0.219;
x0 = 0.5;

Hp = 3;
Q = 1;
S = 1;
R = 1;

A_cal = [A; A^2; A^3]
B_cal = [B zeros(1, 1) zeros(1, 1); A*B B zeros(1, 1); ...
    A^2*B A*B B]
Q_cal = blkdiag(Q, Q, Q);
R_cal = blkdiag(R, R, R);

H = 2*(B_cal'*Q_cal*B_cal + R_cal);
F = 2*A_cal'*Q_cal*B_cal;

u_max = 0.1;
u_min = -0.1;

G = [eye(3); -eye(3)]
h = [u_max*ones(3, 1); -u_min*ones(3, 1)]

U = quadprog(H, x0'*F, G, h)
X = A_cal * x0 + B_cal * U