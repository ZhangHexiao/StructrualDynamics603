%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Newmark-beta Method (Average Acceleration) for linear SDF systems for
% Example with  Earthquake excitation including free vibration
% Course: CIVE 603 - Structural Dynamics - Winter 2018
% Sarven Akcelyan, PhD, McGill University
% Created: Feb 1, 2018
%
% 
%
%% Clear variables, close all figures/tables and clear the command line
clear all; close all; clc;



%% Define Input Parameters of SDF System
Tn = 2.0;                        % Period of SDF system [sec]
wn = 2*pi/Tn;                    % circular frequency (calculated) [rads/sec]
m = 1.0;                         % Seismic mass of SDF system [kips-s^2/in]
k = m/(Tn/(2*pi))^2;             % lateral stiffness (calculated) [kips/in]
z = 0.02;                        % damping ratio of SDF system
c = z*(2*m*wn);                  % damping coeffiecent  (calculated)
u0 = 0;                          % Initial Displacement
v0 = 0;                          % Initial velocity
Tfree= 5*Tn;                     % Free vibration [sec]
g = 386.22;                      % acceleration of gravity [in/sec^2]

%% Load Ground Motion
% Example: El Centro
ElCentro = load('ElCentro.th');  % Load ground motion file
dt_ElCentro = 0.02;              % Ground motion sampling rate
dt = dt_ElCentro;
Ag=[0;ElCentro;zeros(Tfree/dt,1)];
p = -m*Ag*g;                     % External force vector
%% Solver - Newmark Method (Average Acceleration Method)
beta = 1/4; gamma = 1/2;

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
    
    u(i+1) = u(i) + du(i);
    v(i+1) = v(i) + dv(i);
    acc(i+1) = acc(i) + da(i);
    % Time Vector
    Time(i+1) = Time(i) + dt;
end

% Absolute acceleration history = relative acceleration + ground acceleration
ab_acc = acc + Ag' * g;
%% Plot the response of SDF system

% to control the font size in figures
F_SIZE = 12;

% subplots
figure('position',[200 0.0 750 750],'color','white'); 
subplot(4,1,1)
plot(Time,p)
ylabel('a_g(t) [in/sec^2]')
xlabel('Time, t [sec]')
title('El Centro','fontsize',F_SIZE+2)
grid on;
box on

subplot(4,1,2)
plot(Time, u)
ylabel('u(t) [in]')
xlabel('Time, t [sec]')
grid on;
box on

subplot(4,1,3)
plot(Time, v)
ylabel('v(t) [in/sec]')
xlabel('Time, t [sec]')
grid on;
box on

subplot(4,1,4)
plot(Time, ab_acc)
ylabel('a^t(t) [in/sec^{2}]')
xlabel('Time, t [sec]')
grid on;
box on

% Single Figure
figure('position',[200 50 750 300],'color','white');
axes('fontsize',F_SIZE,'fontweight','b');

plot(Time, acc,'-b')
hold on 
plot(Time, ab_acc,'-r','Linewidth',1.5)
xlabel('Time, t [sec]')
ylabel('Acceleration [in/sec^{2}]')
h_legend=legend('Relative a(t)','Absolute a^t(t)')
set(h_legend,'FontSize',10);
grid off;
box on;
axis([0,45,-200,200])  %Limits for axis [X1,X2,Y1,Y2]







