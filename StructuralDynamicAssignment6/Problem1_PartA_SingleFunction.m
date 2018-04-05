%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Assignment-3
% Course: CIVE 603 - Structural Dynamics - Winter 2018
% Hao Shi      260782588
% Hexiao Zhang,266784352
% 
%
%% Clear variables, close all figures/tables and clear the command line
clear all; close all; clc;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Fi = [0.6375,0.9827,1.5778;1.2750,0.9829,-1.1270;1.9125,-1.9642,0.4508];
m_ = 100/386*[1,0,0;0,1,0;0,0,0.5];

%% Define Input Parameters of SDF System
% units [kip,inch,sec]

wn = 12.01;                        
Fi1 = Fi(:,1);
Gama = (Fi1'*m_*[1;1;1])/(Fi1'*m_*Fi1);
L = m_*Fi1;
h = 12*12;
% wn = 25.47;
% Fi1 = Fi(:,2);
% Gama = (Fi2'*m_*[1;1;1])/(Fi2'*m_*Fi2);
% L = m_*Fi1;

% wn = 38.90;
% Fi1 = Fi(:,3);
% Gama = (Fi3'*m_*[1;1;1])/(Fi3'*m_*Fi3);
% L = m_*Fi1;

Tn = 2*pi/wn;                    % circular frequency (calculated) [rads/sec]
m = 1;                           % Seismic mass of SDF system [kips-s^2/in]
k = m*wn^2;             % lateral stiffness (calculated) [kips/in]


z = 0.05;                        % 5% damping ratio of SDF system
c = z*(2*m*wn);                  % damping coefficient (calculated)
u0 = 0;                          % Initial Displacement [in]
v0 = 0;                          % Initial velocity [in/sec]                   
g = 386.22;                      % acceleration of gravity [in/sec^2]

%%%
%% Load Ground Motion%%%%%%
%% Set the type of earthquake 
%%%%% Canoga_Pake %%%%%%%
Canoga = load('Canoga_Park.th');  % Load ground motion file
dt_Canoga = 0.01;                 % Ground motion sampling rate
dt = dt_Canoga;
Ag=[0;Canoga];
PeakOfGM = max(abs(Ag))*g;
p = -m*Ag*g;                      % External force vector      


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

u(i) = phat(i-1)/khat;                   %%Displacement response%%%%%
phat(i) = p(i)-alpha*u(i-1)-beta*u(i);
u(i+1) = phat(i)/khat;
v(i) = (u(i+1)-u(i-1))/(2*dt);
acc(i) = (u(i+1)-2*u(i)+u(i-1))/(dt^2);    %%Acceleration response%%%%%%
u5(i) = Gama*Fi1(3,1)*u(i);                %%Displacement of 5th floor%%
drift3(i) = Gama*(Fi1(3,1)-Fi1(2,1))*u(i); %%Drift of the third floor%%%
drift2(i) = Gama*(Fi1(2,1)-Fi1(1,1))*u(i); %%Drift of the second floor%%%
fshear3(i) = Gama*L(3,1)*Fi1(3,1)*wn^2*u(i);
fshear2(i) = Gama*L(2,1)*Fi1(2,1)*wn^2*u(i);
fshear1(i) = Gama*L(1,1)*Fi1(1,1)*wn^2*u(i);
Mtotal = fshear3(i)*3*h + fshear2(i)*2*h + fshear1(i)*1*h;

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
figure
% to control the font size in figures
F_SIZE = 10;
% Ploth the deformation response
subplot(2,1,1)
plot(Time, u)
ylabel('Dn [in]')
xlabel('Time, t [sec]')
title('Deformation response for Canoga Park','fontsize',F_SIZE+2)
axis([0,20,-3,3])
grid on;
box on

subplot(2,1,2)
plot(Time, acc)
ylabel('An [in/sec^{2}]')
xlabel('Time, t [sec]')
title('Acceleration response for Canoga Park','fontsize',F_SIZE+2)
axis([0,20,-400,400])
grid on;
box on






