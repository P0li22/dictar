clear all
close all
clc

s = tf('s');
Ts = 5e-3;
Gcont = (28.4*s + 119.7) / (s*(s^2+7.15*s+14.7));

zeta = abs(log(0.12)) / sqrt(pi^2 + (log(0.12))^2);
zeta = 1.1 * zeta;
wn = 4.6 / (zeta * 0.15);

%discretization
G = zpk(minreal(c2d(Gcont, Ts, 'zoh')));
[zg, pg, kg] = zpkdata(G, 'v')
A = conv(conv([1 -pg(1)], [1 -pg(2)]), [1 -pg(3)]);
B = kg * conv([1 -zg(1)], [1 -zg(2)]);

%z-p cancellation
%{
pzmap(G);
hold on;
zgrid(zeta, 0, Ts);
%}

A_plus = conv([1 -pg(2)], [1 -pg(3)]);
A_minus = 1;
B_plus = [1 -zg(2)];
B_minus = kg * [1 -zg(1)];
A_dioph = conv(conv(conv([1 -1], [1, -1]), [1 -1]), A_minus);
B_dioph = B_minus;

%poles
p1c = -zeta*wn + j*wn*sqrt(1-zeta^2);
p2c = -zeta*wn - j*wn*sqrt(1-zeta^2);
p3c = -5*zeta*wn;
p4c = -5*zeta*wn;

p1d = exp(p1c*Ts);
p2d = exp(p2c*Ts);
p3d = exp(p3c*Ts);
p4d = exp(p4c*Ts);

Am = poly([p1d, p2d, p3d, p4d]);

M_S = [ [A_dioph(:); 0], [0; A_dioph(:)], ...
    [0; B_dioph(:); 0; 0], [0; 0; B_dioph(:); 0], [0; 0; 0; B_dioph(:)] ]
gamma = [Am'];
theta = M_S \ gamma;
R1 = theta(1:2)';
S1 = theta(3:5)';
S = conv(S1, A_plus);
R = conv(conv(conv([1 -1], [1 -1]), B_plus), R1);
C = zpk(minreal(tf(S, R, Ts)))

% F

W_prime = zpk(minreal((C * G) / (1 + C * G)));
kW = 1;
T_tilde = conv([1 -p3d], [1 -p4d]);
kT = kW * polyval(S1, 1) / (polyval(T_tilde, 1) * dcgain(W_prime));
T1 = kT * T_tilde;

F = zpk(minreal(tf(T1, S1, Ts)))

%performance evaluation

out_ts = sim("settling_time.slx");
figure(1)
plot(out_ts.y.time, out_ts.y.data);
grid on, zoom on, hold on,
xline(0.15, ':r')
yline(1.01, ':r')
yline(0.99, ':r')

out_u = sim("maxU.slx")
figure(2)
plot(out_u.u.time, out_u.u.data)
grid on, zoom on, hold on,
yline(2.5, ':r')









