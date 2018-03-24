function [Time, u, v, acc, ab_acc]=Newmark_AvAcc_Solver(dtG,dtA,AGM,g,Tn,z,Tfree)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This function performs Earthquake Response History Analysis for linear SDF systems
% Utilizing Newmark-beta Method (Average Acceleration)
% Subroutine for Elastic_Response_Spectrum.m
% Course: CIVE 603 - Structural Dynamics - Winter 2018
% Sarven Akcelyan, PhD, McGill University
% Created: Feb 1, 2018
%
% Input: 
% dtG= time step of the record [sec]
% dtA= given time step of the analysis  [sec] (ex. dtG)
% AGM = time history of the external force (record)
% g = gravitational acceleration [other units are computed based on the unit of g]
% Tn = natural period of the SDF
% z = damping ratio (ex. 0.05)
% Tfree = duration of the free vibration after forced vibration (ex. 5*Tn sec)
%
% Output:
% Time = Time Vector [sec]
% u = relative deformation response history 
% v = relative velocity response history 
% acc = relative acceleration response history 
% ab_acc = absolute acceleration response history 
%% Parameters of SDF System
wn = 2*pi/Tn;                    % circular frequency (calculated) [rads/sec]
m = 1.0;                         % Seismic mass of SDF system 
k = m/(Tn/(2*pi))^2;             % lateral stiffness (calculated) 
c = z*(2*m*wn);                  % damping coefficent (calculated) 
u0 = 0;                          % Initial Displacement 
v0 = 0;                          % Initial velocity 

%% Accuracy Check-GM, Interpolation in case dt is different than dtG
% Note: Newmark Average Acceleartion method is unconditionally stable
dt=dtG;
while dt > Tn/10 % For accuracy in small periods
    dt = dt/2; 
end
dt=min(dt,dtA); % final time step of the analysis [sec]
AGMN=interp1(0:dtG:dtG*(length(AGM)),[0;AGM],0:dt:dtG*(length(AGM)));

%% External Force vector
Ag=[AGMN';zeros(round(Tfree/dt),1)];
p = -m*Ag*g;                     % External force vector

%% Solver - Newmark Method (Average Acceleration Method)
beta=0.25; gamma=0.5;
% Step 1: Initial Conditions
p0 = p(1);                       % initial value of external force
a0 = (p0 - c*v0 - k*u0)/m;

kh = k + gamma*c/(beta*dt) + m/(beta*dt^2) ; 
a = m/(beta*dt) + gamma*c/beta;
b = m/(2*beta) + dt*(gamma/(2*beta) - 1)*c;

% Initialize Vectors to be used 
acc(1) = a0;
v(1) = v0;
u(1) = u0;
Time(1) = 0.00;


% Step 2 Calculations for time step i
for i = 1:(length(p)-1); 
    Dph(i) = (p(i+1) - p(i)) + a*v(i) + b *acc(i);
    du(i) = Dph(i)/kh;
    dv(i) = gamma*du(i)/(beta*dt) - gamma*v(i)/beta + dt*(1-gamma/(2*beta))*acc(i);
    da(i) = du(i)/(beta*dt^2) - v(i)/(beta*dt) - acc(i)/(2*beta);
    u(i+1) = u(i) + du(i); % relative displacement
    v(i+1) = v(i) + dv(i); % relative velocity
    acc(i+1) = acc(i) + da(i); % relative acceleration
    Time(i+1) = Time(i) + dt; % Time Vector
end

% Absolute acceleration history = relative acceleration + ground acceleration
ab_acc = acc + Ag' * g;
end






