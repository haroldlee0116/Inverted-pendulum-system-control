clc;
clear all;
close all;
%% Set up parameters
br=0.001;
bp=0.00005;
% syms Rm kt km mr r Jr mp Lp l Jp g Jt br bp
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
%Jp = mp*Lp^2/12;
% Equivalent Viscous Damping Coefficient (N-m-s/rad)
%bp = 0.00005; % damping tuned heuristically to match QUBE-Sero 2 response
%bp = 0.00005; % damping tunedDamping Coefficient (N-m-s/rad)
%br = 0.0033; % damping tuned heuristically to match QUBE-Sero 2 response

% Gravity Constant
g = 9.81;

% syms Rm kt km mr r Jr mp Lp l Jp g Jt br bp

%   Copyright 2013 The MathWorks, Inc.
Jt = Jr*Jp - mp^2*r^2*l^2;
% %
% syms Rm kt km mr r Jr mp Lp l Jp g Jt br bp

% State Space Representation
A = [0 0 1 0;
    0 0 0 1;
    0 mp^2*l^2*r*g/Jt  -br*Jp/Jt  mp*l*r*bp/Jt;
    0  -mp*g*l*Jr/Jt   mp*l*r*br/Jt -Jp*bp/Jt];

% original one A = [0 0 1 0;0 0 0 1;0 mp^2*l^2*r*g/Jt  -br*Jp/Jt  mp*l*r*bp/Jt; 0  -mp*g*l*Jr/Jt mp*l*r*br/Jt -Jr*bp/Jt];
%
Bold = [0; 0; Jp/Jt; -mp*l*r/Jt];
% original one B = [0; 0; Jp/Jt; -mp*l*r/Jt];
C = eye(2,4);
D = zeros(2,1);
% 
% Add actuator dynamics
B = -km * Bold / Rm*40;
A(3,1) = A(3,1) - km*km/Rm*Bold(3);
A(4,1) = A(4,1) + km*km/Rm*Bold(4);

%h=0.5
h=0.01;
Tsim=10;
%% identification plot
sim dowd_id

figure(1)
subplot(2,1,1)
plot(0:2*h:Tsim,alphadata.Data,'LineWidth',2)
hold on
plot(0:2*h:Tsim,alphaiden.Data,':','LineWidth',1.5)
hold off
xlabel('Time(s)')
ylabel('Alpha (rad)')
legend('Collected','P_{r} Exp1')
title('Validation rate for \alpha')

subplot(2,1,2)
plot(0:2*h:Tsim,thetadata.Data,'LineWidth',2)
hold on
plot(0:2*h:Tsim,thetaden.Data,':','LineWidth',1.5)
hold off
xlabel('Time(s)')
ylabel('Theta (rad)')
legend('Collected','P_{r} Exp1')
title('Validation rate for \theta')
%% greyest session
alpha=alphadata.Data;
theta=thetadata.Data;
u=input.Data;
y=[alpha,theta];

data = iddata(y, u,2*h, 'Name', 'Identified system');
data.InputName = 'Voltage';
data.InputUnit = 'V';
data.OutputName = {'Alpha (rad)', 'Theta (rad)'};
data.OutputUnit = {'rad', 'rad'};
data.Tstart = 0;
data.TimeUnit = 's';
data1=data;
save data1
odefun = 'DCMotorODE';

gain=43;
%gain=43;

Ts=h;
parameters = {gain};

fcn_type = 'c';

sys = idgrey(odefun,parameters,fcn_type);

sys.Structure.Parameters(1).Minimum =40;

opt=greyestOptions('EnforceStability',true);
sys=greyest(data,sys, opt);

opt = compareOptions('InitialCondition','zero');
compare(data,sys,Inf,opt)
%% Pole-placement controller
T=ctrb(A,B);
rank(T)
% control specification
zeta=0.7;%damping ratio
wn=4; %natural frequency
% Location of dominant poles along real-axis
sigma=zeta*wn;
% Location of dominant poles along img axis
wd=wn*sqrt(1-zeta^2); %damped natural frequency
% Desired poles 
DP= [-sigma+j*wd, -sigma-j*wd, -45, -40];
% Find control gain using MATLAB pole-placement command
K=acker(A,B,DP);
eig(A-B*K)
%% LQR
r=1;
q=[15 0 0 0; 0 4 0 0; 0 0 0.5 0; 0 0 0 0.2];
K=lqr(A, B, q, r);