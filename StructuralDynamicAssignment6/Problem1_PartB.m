%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Assignment-6
% Course: CIVE 603 - Structural Dynamics - Winter 2018
% Hexiao Zhang,266784352
% Hao Shi      260782588
%
% 
%
%% Clear variables, close all figures/tables and clear the command line
clear all; close all; clc;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Fi = [0.6375,0.9827,1.5778;1.2750,0.9829,-1.1270;1.9125,-1.9642,0.4508];
m_ = 100/386*[1,0,0;0,1,0;0,0,0.5];
Fi1 = Fi(:,1);
Fi2 = Fi(:,2);
Fi3 = Fi(:,3);
L1 = m_*Fi1;
L2 = m_*Fi2;
L3 = m_*Fi3;
Gama = [(Fi1'*m_*[1;1;1])/(Fi1'*m_*Fi1),(Fi2'*m_*[1;1;1])/(Fi2'*m_*Fi2),(Fi3'*m_*[1;1;1])/(Fi3'*m_*Fi3)];
S1 = Gama(1) * m_ *Fi1;
S2 = Gama(2) * m_ *Fi2;
S3 = Gama(3) * m_ *Fi3;
m = 1; 
h = 12*12;
g = 386.22;                      % acceleration of gravity [in/sec^2]
z = 0.05;                        % 5% damping ratio of SDF system  
Tfree = 0;


%% Load Ground Motion%%%%%%
%%%%% Canoga_Pake %%%%%%%
AGM = load('Canoga_Park.th');  % Load ground motion file
dt_Canoga = 0.01;                 % Ground motion sampling rate
dtG = dt_Canoga;
dtA = 0.01;

%% Define Input Parameters of SDF System
% units [kip,inch,sec]
 for i=[1,2,3]
  
switch i
    
    case 1
wn = 12.01;                        
Tn = 2*pi/wn;                    % circular frequency (calculated) [rads/sec]
k = m*wn^2;                      % lateral stiffness (calculated) [kips/in]
[Time1, u1, v1, acc1, ab_acc1, fs1]=CentralDiff_Solver(m,Tn,z,dtG,dtA,AGM,g,Tfree);
u3_1 = Gama(1,i)*Fi1(3,1)*u1;                %%Displacement of 5th floor%%
drift3_1 = (Gama(1,i)*(Fi1(3,1)-Fi1(2,1))*u1)/h;
drift2_1 = (Gama(1,i)*(Fi1(2,1)-Fi1(1,1))*u1)/h;
drift1_1 = (Gama(1,i)*(Fi1(1,1)-0)*u1)/h;
fshear3_1 = S1(3)*wn^2*u1;
fshear2_1 = S1(2)*wn^2*u1+fshear3_1;
fshear1_1 = S1(1)*wn^2*u1+fshear2_1;
Mtotal_1 = fshear3_1*3*h + fshear2_1*2*h + fshear1_1*1*h;


    case 2
wn = 25.47;
Tn = 2*pi/wn;                    % circular frequency (calculated) [rads/sec]
k = m*wn^2;                      % lateral stiffness (calculated) [kips/in]
[Time2, u2, v2, acc2, ab_acc2, fs2]=CentralDiff_Solver(m,Tn,z,dtG,dtA,AGM,g,Tfree);
u3_2= Gama(1,i)*Fi2(3,1)*u2;                %%Displacement of 5th floor%%
drift3_2 = (Gama(1,i)*(Fi2(3,1)-Fi2(2,1))*u2)/h;    
drift2_2 = (Gama(1,i)*(Fi2(2,1)-Fi2(1,1))*u2)/h;
drift1_2 = (Gama(1,i)*(Fi2(1,1)-0)*u2)/h;
fshear3_2 = S2(3)*wn^2*u2;
fshear2_2 = S2(2)*wn^2*u2+fshear3_2;
fshear1_2 = S2(1)*wn^2*u2+fshear2_2;
Mtotal_2 = fshear3_2*3*h + fshear2_2*2*h + fshear1_2*1*h;


    case 3
wn = 38.90;
Tn = 2*pi/wn;                    % circular frequency (calculated) [rads/sec]
k = m*wn^2;                      % lateral stiffness (calculated) [kips/in]
[Time3, u3, v3, acc3, ab_acc3, fs3]=CentralDiff_Solver(m,Tn,z,dtG,dtA,AGM,g,Tfree);
u3_3 = Gama(1,i)*Fi3(3,1)*u3;                %%Displacement of 5th floor%%
drift3_3 = (Gama(1,i)*(Fi3(3,1)-Fi3(2,1))*u3)/h;  
drift2_3 = (Gama(1,i)*(Fi3(2,1)-Fi3(1,1))*u3)/h;
drift1_3 = (Gama(1,i)*(Fi3(1,1)-0)*u3)/h;
fshear3_3 = S3(3)*wn^2*u3;
fshear2_3 = S3(2)*wn^2*u3+fshear3_3;
fshear1_3 = S3(1)*wn^2*u3+fshear2_3;
Mtotal_3 = fshear3_3*3*h + fshear2_3*2*h + fshear1_3*1*h;


end
 end
%%%%%%Response History Analysis%%%%%%%%%%%%%%%
u3 = u3_1+u3_2+u3_3;
drift3 = drift3_1+drift3_2+drift3_3;
drift2 = drift2_1+drift2_2+drift2_3;
drift1 = drift1_1+drift1_2+drift1_3;
fshear3 = fshear3_1+fshear3_2+fshear3_3;
fshear2 = fshear2_1+fshear2_2+fshear2_3;
fshear1 = fshear1_1+fshear1_2+fshear1_3;
Mtotal = Mtotal_1+Mtotal_2+Mtotal_3;

%%%%%Load Sap2000 History Analysis Data%%%%%%%%
u3_sap = load('u3Displacement.txt');
u2_sap = load('u2Displacement.txt');
u1_sap = load('u1Displacement.txt');
f3_sap = load('shearForceThirdFloor.txt');
f2_sap = load('shearForceSecondFloor.txt');
f1_sap = load('shearForceFirstFloor.txt');
M_sap  = load('baseMoment.txt');

Time_sap = u3_sap(:,1);
drift3_sap = (u3_sap(:,2)-u2_sap(:,2))/h;
drift2_sap = (u2_sap(:,2)-u1_sap(:,2))/h;
drift1_sap = (u1_sap(:,2)-0)/h;
%%%laod the data and then plot%%%%




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure(1)
F_SIZE =10;
[maxy,idx] = max(abs(u3));
plot(Time1, u3, Time1(idx),u3(idx),'pr')
hold on
plot (Time_sap,u3_sap(:,2),'--r','Linewidth',1.5)
strmax = [' ',num2str(maxy)];
% text(Time1(idx),-maxy,strmax,'Color','red','FontSize',10);
ylabel('Dn [in]')
xlabel('Time, t [sec]')
title('Response of total roof displacement','fontsize',F_SIZE)
axis([0,20,-5,5])
box on
grid on
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure(2)
F_SIZE =10;
[maxy,idx] = max(abs(drift3));
plot(Time1, drift3, Time1(idx),drift3(idx),'pr')
hold on
plot (Time_sap,drift3_sap,'--r','Linewidth',1.5)
strmax = [' ',num2str(maxy)];
% text(Time1(idx),-maxy,strmax,'Color','red','FontSize',10);
ylabel('Ratio')
xlabel('Time, t [sec]')
title('Response of total third-story drift ratio','fontsize',F_SIZE)
axis([0,20,-0.01,0.01])
box on
grid on
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure(3)
F_SIZE =10;
[maxy,idx] = max(abs(drift2));
plot(Time1, drift2, Time1(idx),drift2(idx),'pr')
hold on
plot (Time_sap,drift2_sap,'--r','Linewidth',1.5)
strmax = [' ',num2str(maxy)];
% text(Time1(idx),-maxy,strmax,'Color','red','FontSize',10);
ylabel('Ratio')
xlabel('Time, t [sec]')
title('Response of total secnd-story drift ratio','fontsize',F_SIZE)
axis([0,20,-0.01,0.01])
box on
grid on
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure(4)
F_SIZE =10;
[maxy,idx] = max(abs(drift1));
plot(Time1, drift1, Time1(idx),drift1(idx),'pr')
hold on
plot (Time_sap,drift1_sap,'--r','Linewidth',1.5)
strmax = [' ',num2str(maxy)];
% text(Time1(idx),-maxy,strmax,'Color','red','FontSize',10);
ylabel('Ratio')
xlabel('Time, t [sec]')
title('Response of total first-story drift ratio','fontsize',F_SIZE)
axis([0,20,-0.01,0.01])
box on
grid on

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure(5)
F_SIZE =10;
[maxy,idx] = max(abs(fshear3));
plot(Time1, fshear3, Time1(idx),fshear3(idx),'pr')
hold on
plot (f3_sap(:,1),2*f3_sap(:,2),'--r','Linewidth',1.5)
strmax = [' ',num2str(maxy)];
% text(Time1(idx),-maxy,strmax,'Color','red','FontSize',10);
ylabel('Shear, f [kips]')
xlabel('Time, t [sec]')
title('Response of total third-story shear force','fontsize',F_SIZE)
axis([0,20,-100,100])
box on
grid on
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure(6)
F_SIZE =10;
[maxy,idx] = max(abs(fshear2));
plot(Time1, fshear2, Time1(idx),fshear2(idx),'pr')
hold on
plot (f2_sap(:,1),2*f2_sap(:,2),'--r','Linewidth',1.5)
strmax = [' ',num2str(maxy)];
% text(Time1(idx),-maxy,strmax,'Color','red','FontSize',10);
ylabel('Shear, f [kips]')
xlabel('Time, t [sec]')
title('Response of total second-story sheasr force','fontsize',F_SIZE)
axis([0,20,-150,150])
box on
grid on
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure(7)
F_SIZE =10;
[maxy,idx] = max(abs(fshear1));
plot(Time1, fshear1, Time1(idx),fshear1(idx),'pr')
hold on
plot (f1_sap(:,1),2*f1_sap(:,2),'--r','Linewidth',1.5)
strmax = [' ',num2str(maxy)];
% text(Time1(idx),-maxy,strmax,'Color','red','FontSize',10);
ylabel('Shear, f [kips]')
xlabel('Time, t [sec]')
title('Response of total first-story shear force','fontsize',F_SIZE)
axis([0,20,-300,300])
box on
grid on
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure(8)
F_SIZE =10;
[maxy,idx] = max(abs(Mtotal));
plot(Time1, Mtotal, Time1(idx),Mtotal(idx),'pr')
hold on
plot (M_sap(:,1),-2*M_sap(:,2),'--r','Linewidth',1.5)
strmax = [' ',num2str(maxy)];
% text(Time1(idx),-maxy,strmax,'Color','red','FontSize',10);
ylabel('Moment, M [kips.in]')
xlabel('Time, t [sec]')
title('Response of total overturning moment','fontsize',F_SIZE)
axis([0,20,-120000,120000])
box on
grid on


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%