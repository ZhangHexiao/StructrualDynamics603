%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Assignment-3
% Course: CIVE 603 - Structural Dynamics - Winter 2018
% Hao Shi      260782588
% Hexiao Zhang 266784352
%
%% Clear variables, close all figures/tables and clear the command line
clear all; close all; clc;


%% Define Input Parameters of SDF System
% units [kip,inch,sec]
Tn = 1.0;                        % Period of SDF system [sec]
wn = 2*pi/Tn;                    % circular frequency (calculated) [rads/sec]
m = 1;                           % Seismic mass of SDF system [kips-s^2/in]
k = m/(Tn/(2*pi))^2;             % lateral stiffness (calculated) [kips/in]
z = 0.0;                         % damping ratio of SDF system
c = z*(2*m*wn);                  % damping coefficient (calculated)
u0 = 0;                          % Initial Displacement [in]
v0 = 0;                          % Initial velocity [in/sec]
Tfree= 3*Tn;                     % Free vibration duration after forced vibration [sec]
g = 386.22;                      % acceleration of gravity [in/sec^2]
%% Free and Force Vibration
% Harmonic Try with and without damping                   
dt = 0.01;
w = 2/3*wn;                           % Excitation frequency- Try 0.01, 0.5, 1, 5, 100*wn (use dt=0.001 for 100*wn)
po = -0.5*g*m;
pt = po*sin(w*[0:dt:1*2*pi/w]);      % 1 cycles of harmonic excitation
%p = [pt];
p = [pt zeros(1,Tfree/dt)]; 


%% Numerical Solution-Central Difference Method
khat = m/dt^2+c/(2*dt);
alpha = m/dt^2-c/(2*dt);
beta = k-2*m/dt^2;
a0 = (p(1)-c*v0-k*u0)/m;
u_1 = u0-dt*v0+dt^2/2*a0;

% Initialize Vectors to be used 

Time(1) = 0.00;%% Beginning point
u(1) = u0;%%beginning point
u_analytical(1) = 0;
phat(1) = p(1)-alpha*u_1-beta*u(1);%% phat0

% Step 2 Calculations for time step i
for i = 2:(length(p)-1)  
    
% Time Vector
Time(i) = Time(i-1) + dt;

if Time(i)<1.5
u_analytical(i) = -(9*g/(40*pi*pi))*(sin(4/3*pi*Time(i))-2/3*sin(2*pi*Time(i)));    
else   u_analytical(i) = -0.3*g/(pi*pi)*sin(2*pi*(Time(i)-1.5));
end

u(i) = phat(i-1)/khat;%%u1
phat(i) = p(i)-alpha*u(i-1)-beta*u(i);
    

end


%% Plot the response of SDF system

% to control the font size in figures
F_SIZE = 10;

figure('position',[200 0.0 750 750],'color','white'); 
subplot(4,1,1)
plot(Time,u)
ylabel('u(t) [in/sec^2]')
xlabel('Time, t [sec]')
title('Problem-1 Part A','fontsize',F_SIZE+2)
grid on;
box on

subplot(4,1,2)
plot(Time, u_analytical)
ylabel('u-analytical(t) [in]')
xlabel('Time, t [sec]')
grid on;
box on
