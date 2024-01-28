clear all
close all 
clc

s = tf('s');
Gcont = 0.2/( (1+s/20)*(1+s/5) );
Ts = 10e-3;

%discretization
G = zpk(minreal(c2d(Gcont, Ts, 'zoh')))
[zg, pg, kg] = zpkdata(G, 'v')
A = conv([1 -pg(1)], [1 -pg(2)]);
B = kg * [1 -zg];

% zeta, wn
s_hat = 0.2;
ts = 1.5;
zeta = abs(log(s_hat))/sqrt(pi^2+(log(s_hat))^2);
zeta = 1.1*zeta
wn = 4.6/(zeta*ts)

%z-p cancellation
%{
pzmap(G);
hold on
zgrid(zeta,0,Ts);
%}
A_plus = A;
A_minus = 1;
B_plus = 1;
B_minus = B;
A_dioph = conv([1 -1], A_minus);
B_dioph = B_minus;

%poles
p1c = -zeta*wn + j*wn*sqrt(1-zeta^2);
p2c = -zeta*wn - j*wn*sqrt(1-zeta^2);
p3c = -5*zeta*wn;

p1d = exp(p1c*Ts);
p2d = exp(p2c*Ts);
p3d = exp(p3c*Ts);

Am = poly([p1d, p2d, p3d]);

%Sylvester Matrix

M_S = [ [A_dioph(:); 0; 0], [0; A_dioph(:); 0], [0; 0; A_dioph(:)]...
    [0; B_dioph(:); 0], [0; 0; B_dioph(:)]; ...
    [10e-3, 10e-3, 10e-3, -1.768e-4, -1.768e-4] ];
gamma = [Am'; 0];
theta = M_S \ gamma;

%C def
R1 = theta(1:3)';
S1 = theta(4:5)';
R = conv([1 -1], R1);
S = conv(A_plus, S1);
C = zpk(minreal(tf(S, R, Ts)))

%F
L = zpk(minreal(C*G));
Wprime = zpk(minreal(L/(1+L)));
kW = 1;
T_tilde = [1 -p3d];
kT = kW * polyval(S1, 1) / (polyval(T_tilde, 1) * dcgain(Wprime));
T1 = kT*T_tilde;
F = zpk(minreal(tf(T1, S1, Ts)))

%performance
out = sim("simulink_2023_05.slx");
figure(1)
plot(out.y.time, out.y.data)
hold on, grid on, zoom on
yline(1.01, ':r'); yline(0.99, ':r'); xline(1.5, ':r');

figure(2)
plot(out.u.time, out.u.data)
hold on, grid on, zoom on
yline(12.5, ':r')






















