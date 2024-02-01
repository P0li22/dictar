clear all 
close all
clc

Ts = 5e-3;
s = tf('s');
Gcont = 10*(s+25) / (s^2+25);

%discretization
G = zpk(c2d(Gcont, Ts, 'zoh'))
[zg, pg, kg] = zpkdata(G, 'v')
A = conv([1 -pg(1)], [1 -pg(2)]);
B = kg*[1 -zg];

%zeta and wn
s_hat = 0.2;
zeta = abs(log(s_hat)) / sqrt(pi^2+(log(s_hat))^2);
wn = 4.6/(zeta*0.3);

%z-p cancS

A_plus = 1;
A_minus = A;
B_plus = [1 -zg];
B_minus = kg;
A_dioph = conv([1 -1], A_minus);
B_dioph = B_minus;

p1c = -zeta*wn + j*wn*sqrt(1-zeta^2);
p2c = -zeta*wn - j*wn*sqrt(1 - zeta^2);
p3c = -5*zeta*wn;

p1d = exp(p1c*Ts);
p2d = exp(p2c*Ts);
p3d = exp(p3c*Ts);

Am = poly([p1d, p2d, p3d])

%Sylvester Matrix
M_S = [ [A_dioph(:)], ...
    [0; B_dioph(:); 0; 0], [0; 0; B_dioph(:); 0], [0; 0; 0; B_dioph(:)] ]
gamma = [Am'];
theta = M_S \ gamma;

R1 = theta(1)';
S1 = theta(2:4)';
S = S1;
R = R1 * conv([1 -1], B_plus);
C = zpk(minreal(tf(S, R, Ts)))

L = zpk(minreal(C*G));
Wprime = zpk(minreal(L / (1+L)))

%performance
out = sim("algebraic_model.slx");
figure(1)
subplot(221)
plot(out.y.time, out.y.data)
hold on, zoom on, grid on,
xlabel('t'), ylabel('y')
xline(0.3, ':r');
yline(1.01, ':r'); yline(0.99, ':r');

subplot(222)
plot(out.u.time, out.u.data);
hold on, zoom on, grid on,
xlabel('t'), ylabel('u')
yline(10, ':r')











