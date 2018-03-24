%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Assignment-3
% Course: CIVE 603 - Structural Dynamics - Winter 2018
% Hao Shi      260782588
% Hexiao Zhang,266784352
% 
%
%% Clear variables, close all figures/tables and clear the command line
clear all; close all; clc;


%% Define Input Parameters of SDF System
% units [kip,inch,sec]
Tn = 2.0;                        % Period of SDF system [sec]
wn = 2*pi/Tn;                    % circular frequency (calculated) [rads/sec]
m = 1;                           % Seismic mass of SDF system [kips-s^2/in]
k = m/(Tn/(2*pi))^2;             % lateral stiffness (calculated) [kips/in]
z = 0.05;                        % 5% damping ratio of SDF system
% z = 0.0;                       % 0% damping ratio of SDF system
c = z*(2*m*wn);                  % damping coefficient (calculated)
u0 = 0;                          % Initial Displacement [in]
v0 = 0;                          % Initial velocity [in/sec]                   
g = 386.22;                      % acceleration of gravity [in/sec^2]

%%%
%% Load Ground Motion%%%%%%
%% Set the type of earthquake

%earthquake = 'EI';

earthquake = 'Canoga'  

%%%% El Centro %%%%%%%%%
switch earthquake
    case 'EI'     
ElCentro = load('ElCentro.th');   % Load ground motion file
dt_ElCentro = 0.02;               % Ground motion sampling rate
dt = dt_ElCentro;
Ag=[0;ElCentro];
PeakOfGM = max(abs(Ag))*g
p = -m*Ag*g;                      % External force vector
%%%%% Canoga_Pake %%%%%%%
    case 'Canoga'
Canoga = load('Canoga_Park.th');  % Load ground motion file
dt_Canoga = 0.01;                 % Ground motion sampling rate
dt = dt_Canoga;
Ag=[0;Canoga];
PeakOfGM = max(abs(Ag))*g;
p = -m*Ag*g;                      % External force vector      
end

%% Numerical Solution-Central Difference Method
khat = m/dt^2+c/(2*dt);
alpha = m/dt^2-c/(2*dt);
beta = k-2*m/dt^2;
acc0 = (p(1)-c*v0-k*u0)/m;
u_1 = u0-dt*v0+dt^2/2*acc0;

%% Step 1 Initialize Vectors to be used 

Time(1) = 0.00;%% Beginning point
u(1) = u0;%%beginning point
v(1) = v0;
acc(1) = acc0;
phat(1) = p(1)-alpha*u_1-beta*u(1);%% phat0
Aground(1) = Ag(1);


%% Step 2 Calculations for time step i
for i = 2:(length(p))    

u(i) = phat(i-1)/khat;%%u1
phat(i) = p(i)-alpha*u(i-1)-beta*u(i);
u(i+1) = phat(i)/khat;
v(i) = (u(i+1)-u(i-1))/(2*dt);
acc(i) = (u(i+1)-2*u(i)+u(i-1))/(dt^2);
% Time Vector
Time(i) = Time(i-1) + dt;
Aground(i) = Ag(i);
if Time(i)>20
   TimeEnd = i;
    break
end       
end

u(:,TimeEnd) = [];

%%absolute acceleration%% 
acc_absolute = acc + Aground * g;

%% Plot the response of SDF system

% to control the font size in figures
F_SIZE = 10;
% Ploth the deformation response
plot(Time, u)
ylabel('u [in/sec^{2}]')
xlabel('Time, t [sec]')
title('Deformation response of 5% damping ration for Canoga Park','fontsize',F_SIZE+2)
%axis([0,20,-20,20])
grid on;
box on








