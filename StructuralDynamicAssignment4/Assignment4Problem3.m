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
m = 1;                           % Seismic mass of SDF system [kips-s^2/in]
z = 0.05;                        % 5% damping ratio of SDF system
u0 = 0;                          % Initial Displacement [in]
v0 = 0;                          % Initial velocity [in/sec]
g = 386.22;                      % acceleration of gravity [in/sec^2]


%%%
%% Load Ground Motion%%%%%%
%%%%% Canoga_Pake %%%%%%%

Canoga = load('Canoga_Park.th');  % Load ground motion file
dt_Canoga = 0.01;                 % Ground motion sampling rate
AGM=[0;Canoga];
dtG =0.01;
dt = 0.0001;
Ag=interp1(0:dtG:dtG*(length(AGM)),[0;AGM],0:dt:dtG*(length(AGM)));

PeakOfGM = max(abs(Ag))*g;
p = -m*Ag*g; % External force vector


T(1)=0.0;
dT=0.05;

% for n = 2:600
  for n = 2:120  
    T(n) = T(n-1) + dT;
    Tn = T(n);                       % Period of SDF system [sec]
    wn = 2*pi/Tn;                    % circular frequency (calculated) [rads/sec]
    c = z*(2*m*wn);                  % damping coefficient (calculated)
    k = m*wn^2;             % lateral stiffness (calculated) [kips/in]
    U0(1)=0;
    %%%%For the elasto-plastic analysis%%%%
        
        %% Numerical Solution-Central Difference Method
        khat = m/dt^2+c/(2*dt);
        alpha = m/dt^2-c/(2*dt);
        acc0 = (p(1)-c*v0-k*u0)/m;
        u_1 = u0-dt*v0+dt^2/2*acc0;
        
        %% Step 1 Initialize Vectors to be used
        
        Time(1) = 0.00;%% Beginning point
        u(1) = u0;%%beginning point
        v(1) = v0;
        beta = k-2*m/(dt^2);
        acc(1) = acc0;
        phat(1) = p(1)-alpha*u_1-beta*u(1);%% phat0
        Aground(1) = Ag(1);   
        for i = 2:(length(p))
            %%%%%%%%%%%%%%%%%%%%%%%%%%
                beta = k-2*m/(dt^2);%%change
                u(i) = phat(i-1)/khat;%%u1
                phat(i) = p(i)-alpha*u(i-1)-beta*u(i);%%change
                u(i+1) = phat(i)/khat;
                v(i) = (u(i+1)-u(i-1))/(2*dt);
                acc(i) = (u(i+1)-2*u(i)+u(i-1))/(dt^2);
                % Time Vector
                Time(i) = Time(i-1) + dt;
                Aground(i) = Ag(i);
                TimeEnd = i;
        end
%       u(:,TimeEnd) = [];
        U0(n) = max(abs(u));
        
        
        Ry(1)=0.5;
        dRy = 0.1;    
    for j = 2:100      
        Ry(j)=Ry(j-1)+dRy;
%       Ry=1;%%%%%%%%%%%%%%%%%%%%%%%
        f0=k*U0(n);
        fy=f0/Ry(j);%%%%%%%%%%%%%%%%%%%%%%%
        uy= U0(n)/Ry(j);%%%%%%%%%%%%%%%%%%%%%%%
        
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
                u(i+1) = phat(i)/khat;
                v(i) = (u(i+1)-u(i-1))/(2*dt);
                acc(i) = (u(i+1)-2*u(i)+u(i-1))/(dt^2);
                f(i)=k*(u(i)-xmin-uy);
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
                f(i)=k*u(i)-k*(xmax-uy);
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
            acc_ab(i) = acc(i) + Aground(i) * g;
        end 
        umax(n)=max(abs(u));
        miu(n) = umax(n)/uy;
        
        if  (miu(n)<(1+0.1))&&(miu(n)>(1-0.1))
            aspec1(n) =max(abs(acc_ab));
            Ry(j);
            Ryy(n) = Ry(j);      
        end
         
        if  (miu(n)<(1.5+0.1))&&(miu(n)>(1.5-0.1))
            aspec1_5(n) =max(abs(acc_ab));
            Ry(j);
            Ryy(n) = Ry(j);
        end       
        
        if  (miu(n)<(2+0.1))&&(miu(n)>(2-0.1))
            aspec2(n) =max(abs(acc_ab));
            Ry(j);
            Ryy(n) = Ry(j);
            
        end        
        
        if  (miu(n)<(4+0.1))&&(miu(n)>(4-0.1))
            aspec4(n) =max(abs(acc_ab));
            Ry(j);
            Ryy(n) = Ry(j);
        end       
        
        if  (miu(n)<(8+0.1))&&(miu(n)>(8-0.1))
            aspec8(n) =max(abs(acc_ab));
            Ry(j);
            Ryy(n) = Ry(j);
        end        
        
        
    end
 end
% aspec(:,TimeEnd) = [];


%% Plot the response of SDF system
figure
% to control the font size in figures
F_SIZE = 10;
% Ploth the deformation response

plot(T, aspec1)
hold on
plot(T, aspec1_5,'Linewidth',2)
hold on
plot(T, aspec2,'Linewidth',2.5)
hold on
plot(T, aspec4,'Linewidth',3)
hold on
plot(T, aspec8,'Linewidth',4)
hold on

ylabel('Absolute Acceleration [in/sec^{2}]')
xlabel('Period, t [sec]')
title('Constant ductility absolute acceleration response spectra','fontsize',F_SIZE+2)
h_legend=legend('u=1','u=1.5','u=2','u=4','u=8');
axis([0.3,12,-inf,inf])
grid on;
box on

% figure
% plot(u,f)
% ylabel('Force [kip]')
% xlabel('Displacement,[inch]')
% %axis([-15,5,-25,35])
% grid on;
% box on







