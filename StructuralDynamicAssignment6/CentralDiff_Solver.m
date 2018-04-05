function [Time, u, v, acc, ab_acc, fs]=CentralDiff_Solver(m,Tn,z,dtG,dtA,AGM,g,Tfree)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This function performs Earthquake Response History Analysis for linear SDF systems
% By utilizing Central Difference Method
% Subroutine for Elastic_Response_Spectrum.m
% Course: CIVE 603 - Structural Dynamics - Winter 2017
% Instructor: Sarven Akcelyan, McGill University
% Created: Feb 19, 2017
%
% Input:
% m = mass of the SDF
% Tn = natural period of the SDF
% z = damping ratio (ex. 0.05)
% dtG= time step of the record [sec]
% dtA= given time step of the analysis  [sec] (ex. dtG)
% AGM = time history of the external force (record)
% g = gravitational acceleration [other units are computed based on the unit of g]
% Tfree = duration of the free vibration after forced vibration (ex. 5*Tn sec)
%
% Output:
% Time = Time Vector [sec]
% u = relative deformation response history
% v = relative velocity response history
% acc = relative acceleration response history
% ab_acc = absolute acceleration response history
% fs = elastic resisting force history
%% Parameters of SDF System
wn = 2*pi/Tn;                    % circular frequency (calculated) [rads/sec]
k = m/(Tn/(2*pi))^2;             % lateral stiffness (calculated)
c = z*(2*m*wn);                  % damping coefficient (calculated)
u0 = 0;                          % Initial Displacement
v0 = 0;                          % Initial velocity

%% Accuracy Check-GM, Interpolation in case dt is different than dtG
% Note: Central difference method is conditionally stable for dt < Tn/pi
dt=dtG;
while dt > Tn/10 % For accuracy in small periods (this is smaller than condition)
    dt = dt/2;
end
dt=min(dt,dtA); % final time step of the analysis [sec]

if dt<dtG
    AGMN=interp1(0:dtG:dtG*(length(AGM)),[0;AGM],0:dt:dtG*(length(AGM)));
else
    AGMN=[0;AGM]';
end
%% External Force vector
 Ag=[AGMN';zeros(round(Tfree/dt),1)];
p = -m*Ag*g;                     % External force vector

%% Solver - Central Difference Method)
% Step 1: Initial Calculations
p0 = p(1);                       % initial value of external force
a0 = (p0 - c*v0 - k*u0)/m;
u_1 = u0 - dt*v0 + dt^2/2*a0;
kh = m/dt^2 + c/(2*dt);
a = m/dt^2 - c/(2*dt);
b = k - 2*m/dt^2;

% Initialize Vectors to be used
Time(1) = 0.00;
acc(1) = a0;
v(1) = v0;
u(1) = u0;
ph(1) = p0 - a*u_1 - b*u0;
u(2) = ph(1)/kh;

% Step 2 Calculations for time step i
for i = 2:(length(p));
    ph(i) = p(i) - a*u(i-1) - b*u(i);
    u(i+1) = ph(i)/kh;
    v(i) = (u(i+1) - u(i-1))/(2*dt);
    acc(i) = (u(i+1) - 2*u(i) + u(i-1))/dt^2;
    Time(i) = Time(i-1) + dt;
end
u(end)=[];  % remove last extra u due to u(i+1);
% Absolute acceleration history = relative acceleration + ground acceleration
ab_acc = acc + Ag' * g;
fs = k*u; % fs = elastic resisting force history
fd = c*v; % fd = damping force history
end






