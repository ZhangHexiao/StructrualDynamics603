%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Assignment-4
% Course: CIVE 603 - Structural Dynamics - Winter 2018
% Hao Shi      260782588
% Hexiao Zhang 266784352
%
%
%% Clear variables, close all figures/tables and clear the command line
clear all; close all; clc;


%% Define Input Parameters of SDF System
% units [kip,inch,sec]
w=200;
h=144;
m = 200/386.22;                  % Seismic mass of SDF system [kips-s^2/in]
Tn = 1.0;                        % Period of SDF system [sec]
Ry = 4;

%Ry = 2;
%Tn = 0.5;                        
%Tn = 1.0;                        
%Tn = 1.0;                        

wn = 2*pi/Tn;                    % circular frequency (calculated) [rads/sec]
k = (wn^2)*m;                       % lateral stiffness (calculated) [kips/in]
z = 0.05;                        % 5% damping ratio of SDF system
c = z*(2*m*wn);                  % damping coefficient (calculated)
u0 = 0;                          % Initial Displacement [in]
v0 = 0;                          % Initial velocity [in/sec]
g = 386.22;                      % acceleration of gravity [in/sec^2]

%%%
%% Load Ground Motion%%%%%%
%% Set the type of earthquake

%%%%% Canoga_Pake %%%%%%%
LomaPrieta = load('LomaPrieta.th');  % Load ground motion file
dt_LomaPrieta = 0.005;               % Ground motion sampling rate
AGM=[0;LomaPrieta];
dtG =0.005;
dt = 0.0001;
Ag=interp1(0:dtG:dtG*(length(AGM)),[0;AGM],0:dt:dtG*(length(AGM)));
%%%%%%%%%%
PeakOfGM = max(abs(Ag))*g;
p = -m*Ag*g; % External force vector
%%%%For the elasto-plastic analysis%%%%

% U0=2.093;     %%By setting Ry=1, we can find that u0=2.09%
U0=2.723;     %%By setting Ry=1, we can find that u0=2.09%
uy =U0/Ry;                           
% uy=10000000;
%% Numerical Solution-Central Difference Method
khat = m/dt^2+c/(2*dt);
alpha = m/dt^2-c/(2*dt);
acc0 = (p(1)-c*v0-k*u0)/m;
u_1 = u0-dt*v0+dt^2/2*acc0;

%% Step 1 Initialize Vectors to be used

Time(1) = 0.00;%% Beginning point
u(1) = u0;%%beginning point
v(1) = v0;
acc(1) = acc0;
beta = k-2*m/dt^2;%%change
phat(1) = p(1)-alpha*u_1-beta*u(1);%% phat0
Aground(1) = Ag(1);
f(1)=0;


%% Step 2 Calculations for time step i
stage =1;
xlim = uy;
xmin= -xlim;

for i = 2:(length(p))
    %%%%%%%%%%%%%%%%%%%%%%%%%%
    if(stage ==1)
        beta = k-2*m/dt^2;%%change
        u(i) = phat(i-1)/khat;%%u1
        phat(i) = p(i)+k*xmin+k*uy-alpha*u(i-1)-beta*u(i);%%change
        u(  i+1) = phat(i)/khat;
        v(i) = (u(i+1)-u(i-1))/(2*dt);
        acc(i) = (u(i+1)-2*u(i)+u(i-1))/(dt^2);
        f(i)=k*(u(i-1)-xmin-uy);
        % Time Vector
        Time(i) = Time(i-1) + dt;
        Aground(i) = Ag(i);
        TimeEnd = i;
        
        if u(i)>xlim
            stage = 2;
        end
        
        if u(i)<xmin
             stage = 4;    
        end

    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%
    if(stage==2)
        
        beta = -2*m/dt^2;%%change
        u(i) = phat(i-1)/khat;%%u1
        phat(i) = p(i)-k*uy-alpha*u(i-1)-beta*u(i);%%change
        u(i+1) = phat(i)/khat;
        v(i) = (u(i+1)-u(i-1))/(2*dt);
        acc(i) = (u(i+1)-2*u(i)+u(i-1))/(dt^2);
        f(i)=k*uy;
        % Time Vector
        Time(i) = Time(i-1) + dt;
        Aground(i) = Ag(i);
        TimeEnd = i;
        
        if v(i)<0
            stage = 3;
            xmax = u(i);
            xlim = xmax - 2*uy;
        end

    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%
    if(stage==3)
        
        beta = k-2*m/dt^2;%%change
        u(i) = phat(i-1)/khat;%%u1
        phat(i) = p(i)+k*(xmax-uy)-alpha*u(i-1)-beta*u(i);%%change
        u(i+1) = phat(i)/khat;
        v(i) = (u(i+1)-u(i-1))/(2*dt);
        acc(i) = (u(i+1)-2*u(i)+u(i-1))/(dt^2);
        f(i)=k*u(i-1)-k*(xmax-uy);
        % Time Vector
        Time(i) = Time(i-1) + dt;
        Aground(i) = Ag(i);
        TimeEnd = i;
        
        if u(i)<xlim
            stage = 4;
        end
        if u(i)>xmax
            stage = 2;
        end
        
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%
    if(stage==4)
        
        beta = -2*m/dt^2;%%change
        u(i) = phat(i-1)/khat;%%u1
        phat(i) = p(i)+k*uy-alpha*u(i-1)-beta*u(i);%%change
        u(i+1) = phat(i)/khat;
        v(i) = (u(i+1)-u(i-1))/(2*dt);
        acc(i) = (u(i+1)-2*u(i)+u(i-1))/(dt^2);
        f(i)=-k*uy;
        % Time Vector
        Time(i) = Time(i-1) + dt;
        Aground(i) = Ag(i);
        TimeEnd = i;
        
        if v(i)>0
            stage = 1;
            xlim = u(i)+2*uy;
            xmin = u(i);
        end
        
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%
end

u(:,TimeEnd) =[];

%%absolute acceleration%%
acc_absolute = acc + Aground * g;
%% Plot the response of SDF system
figure
% to control the font size in figures
F_SIZE = 10;
% Ploth the deformation response
subplot(1,2,1);
plot(Time, u)
ylabel('u [in]')
xlabel('Time, t [sec]')
title('Deformation response of 5% damping ration for Loma Prieta','fontsize',F_SIZE+2)
% axis([0,30,-3,4])
grid on;
box on

subplot(1,2,2);
plot(u,f)
ylabel('Force [kip]')
xlabel('Displacement,[inch]')
title('Force-deformation relation','fontsize',F_SIZE+2)
% axis([-15,5,-25,35])
grid on;
box on







