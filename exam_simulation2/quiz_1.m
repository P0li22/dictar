clear all
close all
clc

z = tf('z', 1);
A = [1 1 0; 0 1 0; 0 0 -1];
eig(A)
Xzi = zpk(minreal(z*inv(z*eye(3)-A)))