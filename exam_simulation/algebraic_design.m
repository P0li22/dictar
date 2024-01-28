clear all
close all
clc

s = tf('s')
Gcont = 16*(s+60)/(s^2+160);

%discretization
Ts = 2e-3;
G = zpk(minreal(c2d(Gcont, Ts, 'zoh')))
[zg, pg, kg] = zpkdata(G, 'v')
B = kg * [1 -zg];
A = conv([1 -pg(1)], [1 -pg(2)]);

%zeta and wn
zeta = abs(log(0.19)) / sqrt(pi^2+(log(0.19))^2);
zeta = 1.3*zeta
wn = 4.6 / (zeta * 0.3)
wn = 0.85*wn;

%z-p cancellation
%{
pzmap(G);
hold on
zgrid(zeta, 0, Ts);
%}
A_plus = 1;
A_minus = A;
B_plus = [1 -zg];
B_minus = kg;
A_dioph = conv([1 -1], A_minus);
B_dioph = B_minus;

%poles
p1c = -zeta*wn + j * wn*sqrt(1-zeta^2);
p2c = -zeta*wn - j * wn*sqrt(1-zeta^2);
p3c = -10*zeta*wn;
p4c = -10*zeta*wn;
p1d = exp(p1c*Ts)
p2d = exp(p2c*Ts)
p3d = exp(p3c*Ts)
p4d = exp(p4c*Ts)

Am = poly([p1d p2d p3d p4d])

%Sylvester matrix
M_S = [ [A_dioph(:); 0], [0; A_dioph(:)],...
    [0; B_dioph(:); 0; 0; 0], [0; 0; B_dioph(:); 0; 0], [0; 0; 0; B_dioph(:); 0], [0; 0; 0; 0; B_dioph(:)];...
    [2.26e-4, 2.26e-4, -0.015, -0.015, -0.015, -0.015] ]
gamma = [Am'; 0]
theta = M_S \ gamma

R1 = theta(1:2)';
S1 = theta(3:6)';
S = S1;
R = conv(conv([1 -1], B_plus), R1);
C = zpk(tf(S, R, Ts))

%F
T_tilde = conv([1 -p3d], [1 -p4d]);
Wprime = zpk(minreal((C*G) / (1 + (C*G))))
kw = 1;
kt = (kw * polyval(S1, 1)) / (polyval(T_tilde, 1) * dcgain(Wprime))
T1 = kt * T_tilde;

F = zpk(minreal(tf(T1, S1, Ts)))

%performance
out = sim("algebraic_design_model.slx");
figure(1)
plot(out.y.time, out.y.data)
hold on, grid on, zoom on,
xline(0.3, ':r')
yline(0.8 + 0.8 * 1/100, ':r')
yline(0.8 - 0.8 * 1/100, ':r')

figure(2)
plot(out.u.time, out.u.data)
hold on, grid on, zoom on,
yline(0.26, ':r')



























