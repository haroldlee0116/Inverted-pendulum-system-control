function [A,B,C,D] = DCMotorODE(gain, Ts)
%DCMOTORODE ODE file representing the dynamics of a DC motor parameterized
%by gain G and time constant Tau.
%
%   [A,B,C,D,K,X0] = DCMOTORODE(G,Tau,Ts) returns the state space matrices
%   of the DC-motor with time-constant Tau and static gain G. The sample
%   time is Ts.
%
%   This file returns continuous-time representation if input argument Ts
%   is zero. If Ts>0, a discrete-time representation is returned.
%
% See also IDGREY, GREYEST.
% % Motor
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
% Equivalent Viscous Damping Coefficient (N-m-s/rad)
br = 0.001; % damping tuned heuristically to match QUBE-Sero 2 response

% % Pendulum Link
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
bp = 0.00005; % damping tuned heuristically to match QUBE-Sero 2 response
% Gravity Constant
g = 9.81;


%   Copyright 2013 The MathWorks, Inc.
Jt = Jr*Jp - mp^2*r^2*l^2;
% 
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
B = -km * Bold / Rm*gain;
A(3,1) = A(3,1) - km*km/Rm*Bold(3);
A(4,1) = A(4,1) + km*km/Rm*Bold(4);

%  if Ts>0 % Sample the model with sample time Ts
%     s = expm([[A B]*Ts; zeros(1,5)]);
%     A = s(1:4,1:4);
%     B = s(1:4,5);
%  end
end