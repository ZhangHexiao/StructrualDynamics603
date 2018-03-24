%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Examples with Free Vibration, Harmonic and Pulse Excitation for linear SDF systems
% Utilizing Newmark-beta Method (Average Acceleration) 
% Course: CIVE 603 - Structural Dynamics - Winter 2018
% Sarven Akcelyan, PhD, McGill University
% Created: Feb 1, 2018
%
% 
%
%% Clear variables, close all figures/tables and clear the command line
clear all; close all; clc;


%% Define Input Parameters of SDF System
% units [kip,inch,sec]
Tn = 1.0;                        % Period of SDF system [sec]
wn = 2*pi/Tn;                    % circular frequency (calculated) [rads/sec]
m = 0.2533;                      % Seismic mass of SDF system [kips-s^2/in]
k = m/(Tn/(2*pi))^2;             % lateral stiffness (calculated) [kips/in]
z = 0.05;                        % damping ratio of SDF system
c = z*(2*m*wn);                  % damping coefficient (calculated)
u0 = 0;                          % Initial Displacement [in]
v0 = 0;                          % Initial velocity [in/sec]
Tfree= 2*Tn;                     % Free vibration duration after forced vibration [sec]
g = 386.22;                      % acceleration of gravity [in/sec^2]

%% Free and Force Vibration
Force='Ha'                          %Fr:Free Vibration, Ha:Harmonic, 
                                    %Re: Rectangular Pulse, HS: Half Sine Pulse
%Free Vibration only
if Force=='Fr'                      
dt = 0.01;
po = 0;
p = zeros(1,Tfree/dt); 

% Harmonic Try with and without damping (10%)
elseif Force=='Ha'                  
dt = 0.01;
w = 1.5*wn;                           % Excitation frequency- Try 0.01, 0.5, 1, 5, 100*wn (use dt=0.001 for 100*wn)
po = 20;
pt = po*sin(w*[0:dt:10*2*pi/w]);      % 10 cycles of harmonic excitation
p = [pt];

%Rectangular Pulse
elseif Force=='Re'                  
dt = 0.001;                           % 0.1,0.01,0.001 to check accuracy with td=2*Tn
po = 20;
td = 1.5*Tn;                            % 0.25, 0.5, 1 1.5, 2
pt = po*ones(1,td/dt);
p = [pt zeros(1,Tfree/dt)];  

%Half Sine Pulse Chopra Example 5.4
elseif Force=='HS'                  
dt = 0.1;
td = 0.6*Tn;
po = 10;
pt = po*sin(pi*[0:dt:td]/td);
p = [pt zeros(1,Tfree/dt)];
end

%% Numerical Solution-Newmark Method (Average Acceleration)
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


%% Plot the response of SDF system

% to control the font size in figures
F_SIZE = 10;

%Relative displacement history
figure('position',[200 50 550 450],'color','white');
subplot(2,1,1)
plot(Time, u)
xlabel('Time, t [sec]','FontSize',F_SIZE)
ylabel('u(t) [in]','FontSize',F_SIZE)
grid on;
box on
% Relative displacement history Normalized
subplot(2,1,2)
usto = po/k;
plot(Time, u/usto,'-b','Linewidth',1)
hold on
plot(Time,p/k/usto,'--r','Linewidth',1.5)
xlabel('Time, t [sec]')
ylabel('u(t)/u_{sto}')
h_legend=legend('Dynamic','Static')
set(h_legend,'FontSize',10);
grid off;
box on

%% Writing an Excel File
FileName=Force;
tableHEAD=[{'t_i [sec]'},{'p_i [kip]'},{'acc_i [in/sec2]'},{'v_i [in/sec]'},{'u_i [in]'}];
tableRES=[Time;p;acc;v;u]';
delete([FileName '.xls'])                   %delete the file if exists
xlswrite(FileName,tableHEAD,1,'A1:E1')
xlswrite(FileName,tableRES,1,'A2:E12')      %First 11 points
