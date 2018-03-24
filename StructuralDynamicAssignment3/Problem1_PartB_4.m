%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Assignment-3
% Course: CIVE 603 - Structural Dynamics - Winter 2018
% Hao Shi      260782588
% Hexiao Zhang,266784352
%
%% Clear variables, close all figures/tables and clear the command line
clear all; close all; clc;

%% Set global variables
m = 1;                              % the System has a unit mass
g = 386.22;                         % acceleration of gravity [in/sec^2]
z = 0.05;                           % 5% damping ratio of SDF system  


%% Load Ground Motion%%%%%%
%% set type of earthquake
earthquake = 'EI';
%earthquake = 'Canoga'; 

%%%% El Centro %%%%%%%%%
switch earthquake
    
    case 'EI'     
ElCentro = load('ElCentro.th');  % Load ground motion file
dt_ElCentro = 0.02;              % Ground motion sampling rate
dt = dt_ElCentro;
Ag=[0;ElCentro];
PeakOfGM = max(abs(Ag))*g
p = -m*Ag*g;                     % External force vector
%%%%% Canoga_Pake %%%%%%%
    case 'Canoga'
Canoga = load('Canoga_Park.th');  % Load ground motion file
dt_Canoga = 0.01;              % Ground motion sampling rate
dt = dt_Canoga;
Ag=[0;Canoga];
PeakOfGM = max(abs(Ag))*g;
p = -m*Ag*g;                     % External force vector      
end

%Get the max response of every natural period
T(1)=0.00;
dT=0.001;


for n = 2:3000
    
T(n) = T(n-1) + dT;
  
%% Define Input Parameters of SDF System
% units [kip,inch,sec]
Tn = T(n);                        % Period of SDF system [sec]
wn = 2*pi/Tn;                       % circular frequency (calculated) [rads/sec]                           % Seismic mass of SDF system [kips-s^2/in]
k = m/(Tn/(2*pi))^2;                % lateral stiffness (calculated) [kips/in]                
c = z*(2*m*wn);                     % damping coefficient (calculated)
u0 = 0;                             % Initial Displacement [in]
v0 = 0;                             % Initial velocity [in/sec]
Tfree= 10*Tn;                          % Free vibration duration after forced vibration [sec]
%%%

%% Numerical Solution-Central Difference Method
khat = m/dt^2+c/(2*dt);
alpha = m/dt^2-c/(2*dt);
beta = k-2*m/dt^2;
acc0 = (p(1)-c*v0-k*u0)/m;
u_1 = u0-dt*v0+dt^2/2*acc0;

% Step 1 Initialize Vectors to be used 

Time(1) = 0.00;%% Beginning point
u(1) = u0;%%beginning point
v(1) = v0;
acc(1) = acc0;
phat(1) = p(1)-alpha*u_1-beta*u(1);%% phat0
Aground(1) = Ag(1);

% Step 2 Calculations for time step i
for i = 2:(length(p)-1)    

u(i) = phat(i-1)/khat;
phat(i) = p(i)-alpha*u(i-1)-beta*u(i);
u(i+1) = phat(i)/khat;
v(i) = (u(i+1)-u(i-1))/(2*dt);
acc(i) = (u(i+1)-2*u(i)+u(i-1))/(dt^2);
Aground(i) = Ag(i);
aspec_ab(i) =acc(i)+Aground(i) * g;

% Time Vector
Time(i) = Time(i-1) + dt;

if Time(i)>20
   TimeEnd = i;
    break
end       
end
u(:,TimeEnd) = [];%%delete the u(i+1) to match the time vector


uspec(n) = max(abs(u));
vspec(n) = max(abs(v));
aspec(n) =max(abs(aspec_ab));
pseudo_v(n) = max(abs(u))*2*pi/Tn;
pseudo_a(n) = max(abs(u))*(2*pi/Tn)^2;


end
pseudo_a(1) = 0;
pseudo_v(1) = 0;
uspec(1) = 0;
vspec(1) = 0;
aspec(1) = PeakOfGM;
%% Plot the response of SDF system
%%%%%%%%%%%%%Canoga%%%%%%%%%%%%%%
switch earthquake
    case 'Canoga'

F_SIZE = 10;
figure('position',[200 50 550 450],'color','white');

subplot(2,1,1)
plot(T, pseudo_a,'-b')
hold on 
plot(T, aspec,'-r','Linewidth',2)
xlabel('Natural period, Tn [sec]','FontSize',F_SIZE)
ylabel('Acceleration [in/sec^{2}]','FontSize',F_SIZE)
axis([0.1,3,0,600])
title('Absolute and pseudo acceleration of Canoga','fontsize',F_SIZE+2)
h_legend=legend('Absolute acceleration','pseudo acceleration');
grid on;
box on


subplot(2,1,2)
plot(T, pseudo_v,'-b')
hold on 
plot(T, vspec,'-r','Linewidth',2)
xlabel('Natural period, Tn [sec]','FontSize',F_SIZE)
ylabel('Velocity [in/sec]','FontSize',F_SIZE)
axis([0.1,3,0,60]);
title('Relative and pseudo velocity of Canoga','fontsize',F_SIZE+2)
h_legend=legend('Relative velocity','pseudo velocity');
grid on;
box on
%%%%%%%%%%%EI%%%%%%%%%%%%%%%
    case 'EI'
F_SIZE = 10;
figure('position',[200 50 550 450],'color','white');

subplot(2,1,1)
plot(T, pseudo_a,'-b')
hold on 
plot(T, aspec,'-r','Linewidth',2)
xlabel('Natural period, Tn [sec]','FontSize',F_SIZE)
ylabel('Acceleration [in/sec^{2}]','FontSize',F_SIZE)
axis([0.1,3,0,400])
title('Absolute and pseudo acceleration of EI','fontsize',F_SIZE+2)
h_legend=legend('Absolute acceleration','pseudo acceleration');
grid on;
box on


subplot(2,1,2)
plot(T, pseudo_v,'-b')
hold on 
plot(T, vspec,'-r','Linewidth',2)
xlabel('Natural period, Tn [sec]','FontSize',F_SIZE)
ylabel('Velocity [in/sec]','FontSize',F_SIZE)
axis([0.1,3,0,40]);
title('Relative and pseudo velocity of EI','fontsize',F_SIZE+2)
h_legend=legend('Relative velocity','pseudo velocity');
grid on;
box on

end

