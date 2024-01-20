clear all
close all
clc

s = tf('s');
Gcont = 1 / ( s * (1+s) );
Ts = 0.15;
G = zpk( c2d(Gcont, Ts, 'zoh') )
%zeta = 0.591;
zeta = 0.591 * 1.15;
wn = 1.208;
%{
pzmap(G);
hold on;
zgrid(zeta, 0, Ts);
%}

%poles
p1c = -zeta*wn + j*wn*sqrt(1-zeta^2);
p2c = -zeta*wn - j*wn*sqrt(1-zeta^2);
p3c = -5*zeta*wn;
p1d = exp(p1c * Ts);
p2d = exp(p2c * Ts);
p3d = exp(p3c * Ts);
Am = poly([p1d, p2d, p3d]);

%z-p cancellation
[zg, pg, kg] = zpkdata(G, 'v')
A = conv([1 -pg(1)], [1 -pg(2)])
B = kg * [1 -zg]
A_plus = [1 -pg(2)]
A_minus = 1
B_plus = 1
B_minus = B

%M_s
A_dioph = conv([1 -1], A_minus)
B_dioph = B_minus

M_s = [ [A_dioph(:); 0; 0], [0; A_dioph(:); 0], [0; 0; A_dioph(:)], ...
    [0; B_dioph(:); 0], [0; 0; B_dioph(:)];...
    [1, 1, 1, -0.119, -0.119] ]
Gamma = [Am'; 0];
theta = M_s \ Gamma

R1 = theta(1:3)';
S1 = theta(4:5)';
C = zpk( tf(conv(A_plus, S1), R1, Ts) )

% performance evaluation
out = sim("x1dof_sim.slx");
figure(1)
title("y")
plot(out.y.Time, out.y.Data)
grid on, hold on
figure(2)
title("u")
plot(out.u.Time, out.u.Data)
grid on, hold on


