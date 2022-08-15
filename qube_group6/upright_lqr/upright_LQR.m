clc;
clear all;
close all;
%% Set up parameters
br=0.001;
bp=0.00005;
% Resistance
Rm = 8.4;
%Rm = 12;
% Current-torque (N-m/A)
kt = 0.042;
% Back-emf constant (V-s/rad)
km = 0.042;

% % Rotary Arm
% Mass (kg)
mr = 0.095;
% Total length (m)
r = 0.085;

% Moment of inertia about pivot (kg-m^2)
Jr = mr*r^2/3;
%Jr = mr*r^2/12;
% Equivalent Viscous % % Pendulum Link
% Mass (kg)
mp = 0.024;
% Total length (m)
Lp = 0.129;
% Pendulum center of mass (m)
l = Lp/2;
% Moment of inertia about pivot (kg-m^2)
Jp = mp*Lp^2/3;
% Equivalent Viscous Damping Coefficient (N-m-s/rad)
% Gravity Constant
g = 9.81;
%   Copyright 2013 The MathWorks, Inc.
Jt = Jr*Jp - mp^2*r^2*l^2;

h=0.01;

% State Space Representation
A = [0 0 1 0;
    0 0 0 1;
    0 mp^2*l^2*r*g/Jt  -br*Jp/Jt  -mp*l*r*bp/Jt;
    0  mp*g*l*Jr/Jt   -mp*l*r*br/Jt -Jp*bp/Jt];

Bold = [0; 0; Jp/Jt; -mp*l*r/Jt];
C = eye(2,4);
D = zeros(2,1);
%% Swing up Control
u_max=10; %Input saturation
ke=0.1; %
%% LQR
gain=45;
B = -km * Bold / Rm*gain;
B(4)=-B(4);
A(3,1) = A(3,1) - km*km/Rm*Bold(3);
A(4,1) = A(4,1) + km*km/Rm*Bold(4);

r_lqr=5000000000000; % update
q_lqr=[5000000000 0 0 0; 0 1000 0 0; 0 0 1000000 0; 0 0 0 100];

K=lqr(A, B, q_lqr, r_lqr)
eig(A-B*K)