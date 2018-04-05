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
u5_1 = Gama(1,i)*Fi1(3,1)*u1;                %%Displacement of 5th floor%%
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
u5_2= Gama(1,i)*Fi2(3,1)*u2;                %%Displacement of 5th floor%%
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
u5_3 = Gama(1,i)*Fi3(3,1)*u3;                %%Displacement of 5th floor%%
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
u5 = u5_1+u5_2+u5_3;
drift3 = drift3_1+drift3_2+drift3_3;
drift2 = drift2_1+drift2_2+drift2_3;
drift1 = drift1_1+drift1_2+drift1_3;
fshear3 = fshear3_1+fshear3_2+fshear3_3;
fshear2 = fshear2_1+fshear2_2+fshear2_3;
fshear1 = fshear1_1+fshear1_2+fshear1_3;
Mtotal = Mtotal_1+Mtotal_2+Mtotal_3;




%%%%%%Response Spectra Analysis%%%%%%%%%%%%%%%%

%%%%ABSSUM%%%%%%%%%%%%
u5_ABS = max(abs(u5_1))+max(abs(u5_2))+max(abs(u5_3));
drift3_ABS = max(abs(drift3_1))+max(abs(drift3_2))+max(abs(drift3_3));
drift2_ABS = max(abs(drift2_1))+max(abs(drift2_2))+max(abs(drift2_3));
drift1_ABS = max(abs(drift1_1))+max(abs(drift1_2))+max(abs(drift1_3));
fshear3_ABS = max(abs(fshear3_1))+max(abs(fshear3_2))+max(abs(fshear3_3));
fshear2_ABS = max(abs(fshear2_1))+max(abs(fshear2_2))+max(abs(fshear2_3));
fshear1_ABS = max(abs(fshear1_1))+max(abs(fshear1_2))+max(abs(fshear1_3));
Mtotal_ABS = max(abs(Mtotal_1))+max(abs(Mtotal_2))+max(abs(Mtotal_3));

%%%%SRSS%%%%%%%%%%%%
u5_SRSS = (max(abs(u5_1))^2+max(abs(u5_2))^2+max(abs(u5_3))^2)^0.5;
drift3_SRSS = (max(abs(drift3_1))^2+max(abs(drift3_2))^2+max(abs(drift3_3))^2)^0.5;
drift2_SRSS = (max(abs(drift2_1))^2+max(abs(drift2_2))^2+max(abs(drift2_3))^2)^0.5;
drift1_SRSS = (max(abs(drift1_1))^2+max(abs(drift1_2))^2+max(abs(drift1_3))^2)^0.5;
fshear3_SRSS = (max(abs(fshear3_1))^2+max(abs(fshear3_2))^2+max(abs(fshear3_3))^2)^0.5;
fshear2_SRSS = (max(abs(fshear2_1))^2+max(abs(fshear2_2))^2+max(abs(fshear2_3))^2)^0.5;
fshear1_SRSS = (max(abs(fshear1_1))^2+max(abs(fshear1_2))^2+max(abs(fshear1_3))^2)^0.5;
Mtotal_SRSS = (max(abs(Mtotal_1))^2+max(abs(Mtotal_2))^2+max(abs(Mtotal_3))^2)^0.5;

%%%%SRSS%%%%%%%%%%%%
w=[12.01,25.47,38.90];

p = zeros(3,3);
for x =1:3
    for y=1:3
        p(x,y)=(0.05^2*(1+w(x)/w(y))^2)/((1-w(x)/w(y))^2+4*0.05^2*w(x)/w(y));
    end
end

u5_CQC =(u5_SRSS^2 + 2*p(1,2)*max(abs(u5_1))*max(abs(u5_2)) + 2*p(1,3)*max(abs(u5_1))*max(abs(u5_3))+ 2*p(2,3)*max(abs(u5_2))*max(abs(u5_3)))^0.5;
drift3_CQC =(drift3_SRSS^2 + 2*p(1,2)*max(abs(drift3_1))*max(abs(drift3_2)) + 2*p(1,3)*max(abs(drift3_1))*max(abs(drift3_3))+ 2*p(2,3)*max(abs(drift3_2))*max(abs(drift3_3)))^0.5;
drift2_CQC =(drift2_SRSS^2 + 2*p(1,2)*max(abs(drift2_1))*max(abs(drift2_2)) + 2*p(1,3)*max(abs(drift2_1))*max(abs(drift2_3))+ 2*p(2,3)*max(abs(drift2_2))*max(abs(drift2_3)))^0.5;
drift1_CQC =(drift1_SRSS^2 + 2*p(1,2)*max(abs(drift1_1))*max(abs(drift1_2)) + 2*p(1,3)*max(abs(drift1_1))*max(abs(drift1_3))+ 2*p(2,3)*max(abs(drift1_2))*max(abs(drift1_3)))^0.5;
fshear3_CQC =(fshear3_SRSS^2 + 2*p(1,2)*max(abs(fshear3_1))*max(abs(fshear3_2)) + 2*p(1,3)*max(abs(fshear3_1))*max(abs(fshear3_3))+ 2*p(2,3)*max(abs(fshear3_2))*max(abs(fshear3_3)))^0.5;
fshear2_CQC =(fshear2_SRSS^2 + 2*p(1,2)*max(abs(fshear2_1))*max(abs(fshear2_2)) + 2*p(1,3)*max(abs(fshear2_1))*max(abs(fshear2_3))+ 2*p(2,3)*max(abs(fshear2_2))*max(abs(fshear2_3)))^0.5;
fshear1_CQC =(fshear1_SRSS^2 + 2*p(1,2)*max(abs(fshear1_1))*max(abs(fshear1_2)) + 2*p(1,3)*max(abs(fshear1_1))*max(abs(fshear1_3))+ 2*p(2,3)*max(abs(fshear1_2))*max(abs(fshear1_3)))^0.5;
Mtotal_CQC =(Mtotal_SRSS^2 + 2*p(1,2)*max(abs(Mtotal_1))*max(abs(Mtotal_2)) + 2*p(1,3)*max(abs(Mtotal_1))*max(abs(Mtotal_3))+ 2*p(2,3)*max(abs(Mtotal_2))*max(abs(Mtotal_3)))^0.5;

%%%%%%Result%%%%%%%%%%
format longG
Result = [u5_ABS,drift3_ABS, drift2_ABS,drift1_ABS,fshear3_ABS,fshear2_ABS, fshear1_ABS,Mtotal_ABS;
          u5_SRSS,drift3_SRSS,drift2_SRSS,drift1_SRSS,fshear3_SRSS,fshear2_SRSS,fshear1_SRSS,Mtotal_SRSS;
          u5_CQC,drift3_CQC,drift2_CQC,drift1_CQC,fshear3_CQC,fshear2_CQC,fshear1_CQC,Mtotal_CQC;
          max(abs(u5)),max(abs(drift3)),max(abs(drift2)),max(abs(drift1)),max(abs(fshear3)),max(abs(fshear2)),max(abs(fshear1)),max(abs(Mtotal))]   









%% Plot the response of SDF system
%%%%%%%%%%%%%roof displacement%%%%%%%%%%%%%%
figure(1)
% to control the font size in figures
F_SIZE = 10;
% Ploth the deformation response
subplot(4,1,1)
[maxy,idx] = max(abs(u5_1));
plot(Time1, u5_1, Time1(idx),u5_1(idx),'pr')
strmax = [' ',num2str(maxy)];
text(Time1(idx),-maxy,strmax,'Color','red','FontSize',10);
ylabel('Dn [in]')
xlabel('Time, t [sec]')
title('Response of roof displacement of Mode 1','fontsize',F_SIZE)
axis([0,20,-5,5])


subplot(4,1,2)
[maxy,idx] = max(abs(u5_2));
plot(Time2, u5_2, Time2(idx),u5_2(idx),'pr')
strmax = [' ',num2str(maxy)];
text(Time2(idx),maxy,strmax,'Color','red','FontSize',10);
ylabel('Dn [in]')
xlabel('Time, t [sec]')
title('Response of roof displacement of Mode 2','fontsize',F_SIZE)
axis([0,20,-5,5])


subplot(4,1,3)
[maxy,idx] = max(abs(u5_3));
plot(Time3, u5_3, Time3(idx),u5_3(idx),'pr')
strmax = [' ',num2str(maxy)];
text(Time3(idx),maxy,strmax,'Color','red','FontSize',10);
ylabel('Dn [in]')
xlabel('Time, t [sec]')
title('Response of roof displacement of Mode 3','fontsize',F_SIZE)
axis([0,20,-5,5])


subplot(4,1,4)
[maxy,idx] = max(abs(u5));
plot(Time1, u5, Time1(idx),u5(idx),'pr')
strmax = [' ',num2str(maxy)];
text(Time1(idx),-maxy,strmax,'Color','red','FontSize',10);
ylabel('Dn [in]')
xlabel('Time, t [sec]')
title('Response of total roof displacement','fontsize',F_SIZE)
axis([0,20,-5,5])




%%%%%%%%%%%%%%%%%%%Story Drift Ratio for story 3%%%%%%%%%%%%%%%%%%%
figure(2)
% to control the font size in figures
F_SIZE = 10;
% Ploth the deformation response
subplot(4,1,1)
[maxy,idx] = max(abs(drift3_1));
plot(Time1, drift3_1, Time1(idx),drift3_1(idx),'pr')
strmax = [' ',num2str(maxy)];
text(Time1(idx),-maxy,strmax,'Color','red','FontSize',10);
ylabel('Ratio')
xlabel('Time, t [sec]')
title('Response of third-story drift ratio of Mode 1','fontsize',F_SIZE)
axis([0,20,-0.01,0.01])


subplot(4,1,2)
[maxy,idx] = max(abs(drift3_2));
plot(Time2, drift3_2, Time2(idx),drift3_2(idx),'pr')
strmax = [' ',num2str(maxy)];
text(Time2(idx),-maxy,strmax,'Color','red','FontSize',10);
ylabel('Ratio')
xlabel('Time, t [sec]')
title('Response of third-story drift ratio of Mode 2','fontsize',F_SIZE)
axis([0,20,-0.01,0.01])

subplot(4,1,3)
[maxy,idx] = max(abs(drift3_3));
plot(Time3, drift3_3, Time3(idx),drift3_3(idx),'pr')
strmax = [' ',num2str(maxy)];
text(Time3(idx),maxy,strmax,'Color','red','FontSize',10);
ylabel('Ratio')
xlabel('Time, t [sec]')
title('Response of third-story drift ratio of Mode 3','fontsize',F_SIZE)
axis([0,20,-0.01,0.01])


subplot(4,1,4)
[maxy,idx] = max(abs(drift3));
plot(Time1, drift3, Time1(idx),drift3(idx),'pr')
strmax = [' ',num2str(maxy)];
text(Time1(idx),-maxy,strmax,'Color','red','FontSize',10);
ylabel('Ratio')
xlabel('Time, t [sec]')
title('Response of total third-story drift ratio','fontsize',F_SIZE)
axis([0,20,-0.01,0.01])

%%%%%%%%%%%%%%%%%%%Story Drift Ratio for story 2%%%%%%%%%%%%%%%%%%%
figure(3)
% to control the font size in figures
F_SIZE = 10;
% Ploth the deformation response
subplot(4,1,1)
[maxy,idx] = max(abs(drift2_1));
plot(Time1, drift2_1, Time1(idx),drift2_1(idx),'pr')
strmax = [' ',num2str(maxy)];
text(Time1(idx),-maxy,strmax,'Color','red','FontSize',10);
ylabel('Ratio')
xlabel('Time, t [sec]')
title('Response of second-story drift ratio of Mode 1','fontsize',F_SIZE)
axis([0,20,-0.01,0.01])


subplot(4,1,2)
[maxy,idx] = max(abs(drift2_2));
plot(Time2, drift2_2, Time2(idx),drift2_2(idx),'pr')
strmax = [' ',num2str(maxy)];
text(Time2(idx),maxy,strmax,'Color','red','FontSize',10);
ylabel('Ratio')
xlabel('Time, t [sec]')
title('Response of second-story drift ratio of Mode 2','fontsize',F_SIZE)
axis([0,20,-0.01,0.01])

subplot(4,1,3)
[maxy,idx] = max(abs(drift2_3));
plot(Time3, drift2_3, Time3(idx),drift2_3(idx),'pr')
strmax = [' ',num2str(maxy)];
text(Time3(idx),maxy,strmax,'Color','red','FontSize',10);
ylabel('Ratio')
xlabel('Time, t [sec]')
title('Response of second-story drift ratio of Mode 3','fontsize',F_SIZE)
axis([0,20,-0.01,0.01])


subplot(4,1,4)
[maxy,idx] = max(abs(drift2));
plot(Time1, drift2, Time1(idx),drift2(idx),'pr')
strmax = [' ',num2str(maxy)];
text(Time1(idx),-maxy,strmax,'Color','red','FontSize',10);
ylabel('Ratio')
xlabel('Time, t [sec]')
title('Response of total second-story drift ratio','fontsize',F_SIZE)
axis([0,20,-0.01,0.01])

%%%%%%%%%%%%%%%%%%%Story Drift Ratio for story 1%%%%%%%%%%%%%%%%%%%
figure(4)
% to control the font size in figures
F_SIZE = 10;
% Ploth the deformation response
subplot(4,1,1)
[maxy,idx] = max(abs(drift1_1));
plot(Time1, drift1_1, Time1(idx),drift1_1(idx),'pr')
strmax = [' ',num2str(maxy)];
text(Time1(idx),-maxy,strmax,'Color','red','FontSize',10);
ylabel('Ratio')
xlabel('Time, t [sec]')
title('Response of first-story drift ratio of Mode 1','fontsize',F_SIZE)
axis([0,20,-0.01,0.01])


subplot(4,1,2)
[maxy,idx] = max(abs(drift1_2));
plot(Time2, drift1_2, Time2(idx),drift1_2(idx),'pr')
strmax = [' ',num2str(maxy)];
text(Time2(idx),maxy,strmax,'Color','red','FontSize',10);
ylabel('Ratio')
xlabel('Time, t [sec]')
title('Response of first-story drift ratio of Mode 2','fontsize',F_SIZE)
axis([0,20,-0.01,0.01])

subplot(4,1,3)
[maxy,idx] = max(abs(drift1_3));
plot(Time3, drift1_3, Time3(idx),drift1_3(idx),'pr')
strmax = [' ',num2str(maxy)];
text(Time3(idx),maxy,strmax,'Color','red','FontSize',10);
ylabel('Ratio')
xlabel('Time, t [sec]')
title('Response of first-story drift ratio of Mode 3','fontsize',F_SIZE)
axis([0,20,-0.01,0.01])


subplot(4,1,4)
[maxy,idx] = max(abs(drift1));
plot(Time1, drift1, Time1(idx),drift1(idx),'pr')
strmax = [' ',num2str(maxy)];
text(Time1(idx),-maxy,strmax,'Color','red','FontSize',10);
ylabel('Ratio')
xlabel('Time, t [sec]')
title('Response of total first-story drift ratio','fontsize',F_SIZE)
axis([0,20,-0.01,0.01])

%%%%%%%%%%%%%%%%%%%Story Shear for story 3%%%%%%%%%%%%%%%%%%%
figure(5)
% to control the font size in figures
F_SIZE = 10;
% Ploth the deformation response
subplot(4,1,1)
[maxy,idx] = max(abs(fshear3_1));
plot(Time1, fshear3_1, Time1(idx),fshear3_1(idx),'pr')
strmax = [' ',num2str(maxy)];
text(Time1(idx),-maxy,strmax,'Color','red','FontSize',10);
ylabel('Shear, f [kips]')
xlabel('Time, t [sec]')
title('Response of third-story shear of Mode 1','fontsize',F_SIZE)
axis([0,20,-100,100])


subplot(4,1,2)
[maxy,idx] = max(abs(fshear3_2));
plot(Time2, fshear3_2, Time2(idx),fshear3_2(idx),'pr')
strmax = [' ',num2str(maxy)];
text(Time2(idx),-maxy,strmax,'Color','red','FontSize',10);
ylabel('Shear, f [kips]')
xlabel('Time, t [sec]')
title('Response of third-story shear of Mode 2','fontsize',F_SIZE)
axis([0,20,-100,100])

subplot(4,1,3)
[maxy,idx] = max(abs(fshear3_3));
plot(Time3, fshear3_3, Time3(idx),fshear3_3(idx),'pr')
strmax = [' ',num2str(maxy)];
text(Time3(idx),maxy,strmax,'Color','red','FontSize',10);
ylabel('Shear, f [kips]')
xlabel('Time, t [sec]')
title('Response of third-story shear of Mode 3','fontsize',F_SIZE)
axis([0,20,-100,100])


subplot(4,1,4)
[maxy,idx] = max(abs(fshear3));
plot(Time1, fshear3, Time1(idx),fshear3(idx),'pr')
strmax = [' ',num2str(maxy)];
text(Time1(idx),-maxy,strmax,'Color','red','FontSize',10);
ylabel('Shear, f [kips]')
xlabel('Time, t [sec]')
title('Response of total third-story drift ratio','fontsize',F_SIZE)
axis([0,20,-100,100])

%%%%%%%%%%%%%%%%%%%Story Shear for story 2%%%%%%%%%%%%%%%%%%%
figure(6)
% to control the font size in figures
F_SIZE = 10;
% Ploth the deformation response
subplot(4,1,1)
[maxy,idx] = max(abs(fshear2_1));
plot(Time1, fshear2_1, Time1(idx),fshear2_1(idx),'pr')
strmax = [' ',num2str(maxy)];
text(Time1(idx),-maxy,strmax,'Color','red','FontSize',10);
ylabel('Shear, f [kips]')
xlabel('Time, t [sec]')
title('Response of second-story shear of Mode 1','fontsize',F_SIZE)
axis([0,20,-150,150])


subplot(4,1,2)
[maxy,idx] = max(abs(fshear2_2));
plot(Time2, fshear2_2, Time2(idx),fshear2_2(idx),'pr')
strmax = [' ',num2str(maxy)];
text(Time2(idx),maxy,strmax,'Color','red','FontSize',10);
ylabel('Shear, f [kips]')
xlabel('Time, t [sec]')
title('Response of second-story shear of Mode 2','fontsize',F_SIZE)
axis([0,20,-150,150])

subplot(4,1,3)
[maxy,idx] = max(abs(fshear2_3));
plot(Time3, fshear2_3, Time3(idx),fshear2_3(idx),'pr')
strmax = [' ',num2str(maxy)];
text(Time3(idx),maxy,strmax,'Color','red','FontSize',10);
ylabel('Shear, f [kips]')
xlabel('Time, t [sec]')
title('Response of second-story shear of Mode 3','fontsize',F_SIZE)
axis([0,20,-150,150])


subplot(4,1,4)
[maxy,idx] = max(abs(fshear2));
plot(Time1, fshear2, Time1(idx),fshear2(idx),'pr')
strmax = [' ',num2str(maxy)];
text(Time1(idx),-maxy,strmax,'Color','red','FontSize',10);
ylabel('Shear, f [kips]')
xlabel('Time, t [sec]')
title('Response of total second-story drift ratio','fontsize',F_SIZE)
axis([0,20,-150,150])


%%%%%%%%%%%%%%%%%%%Story Shear for story 1%%%%%%%%%%%%%%%%%%%
figure(7)
% to control the font size in figures
F_SIZE = 10;
% Ploth the deformation response
subplot(4,1,1)
[maxy,idx] = max(abs(fshear1_1));
plot(Time1, fshear1_1, Time1(idx),fshear1_1(idx),'pr')
strmax = [' ',num2str(maxy)];
text(Time1(idx),-maxy,strmax,'Color','red','FontSize',10);
ylabel('Shear, f [kips]')
xlabel('Time, t [sec]')
title('Response of first-story shear of Mode 1','fontsize',F_SIZE)
axis([0,20,-300,300])



subplot(4,1,2)
[maxy,idx] = max(abs(fshear1_2));
plot(Time2, fshear1_2, Time2(idx),fshear1_2(idx),'pr')
strmax = [' ',num2str(maxy)];
text(Time2(idx),-maxy,strmax,'Color','red','FontSize',10);
ylabel('Shear, f [kips]')
xlabel('Time, t [sec]')
title('Response of first-story shear of Mode 2','fontsize',F_SIZE)
axis([0,20,-300,300])


subplot(4,1,3)
[maxy,idx] = max(abs(fshear1_3));
plot(Time3, fshear1_3, Time3(idx),fshear1_3(idx),'pr')
strmax = [' ',num2str(maxy)];
text(Time3(idx),maxy,strmax,'Color','red','FontSize',10);
ylabel('Shear, f [kips]')
xlabel('Time, t [sec]')
title('Response of first-story shear of Mode 3','fontsize',F_SIZE)
axis([0,20,-300,300])



subplot(4,1,4)
[maxy,idx] = max(abs(fshear1));
plot(Time1, fshear1, Time1(idx),fshear1(idx),'pr')
strmax = [' ',num2str(maxy)];
text(Time1(idx),-maxy,strmax,'Color','red','FontSize',10);
ylabel('Shear, f [kips]')
xlabel('Time, t [sec]')
title('Response of total first-story drift ratio','fontsize',F_SIZE)
axis([0,20,-300,300])


%%%%%%%%%%%%%%%%%%%Overturning moment%%%%%%%%%%%%%%%%%%%
figure(8)
% to control the font size in figures
F_SIZE = 10;
% Ploth the deformation response
subplot(4,1,1)
[maxy,idx] = max(abs(Mtotal_1));
plot(Time1, Mtotal_1, Time1(idx),Mtotal_1(idx),'pr')
strmax = [' ',num2str(maxy)];
text(Time1(idx),-maxy,strmax,'Color','red','FontSize',10);
ylabel('Moment, M [kips.in]')
xlabel('Time, t [sec]')
title('Response of overturning moment of Mode 1','fontsize',F_SIZE)
 axis([0,20,-120000,120000])



subplot(4,1,2)
[maxy,idx] = max(abs(Mtotal_2));
plot(Time2, Mtotal_2, Time2(idx),Mtotal_2(idx),'pr')
strmax = [' ',num2str(maxy)];
text(Time2(idx),maxy,strmax,'Color','red','FontSize',10);
ylabel('Moment, M [kips.in]')
xlabel('Time, t [sec]')
title('Response of overturning moment of Mode 2','fontsize',F_SIZE)
 axis([0,20,-120000,120000])


subplot(4,1,3)
[maxy,idx] = max(abs(Mtotal_3));
plot(Time3, Mtotal_3, Time3(idx),Mtotal_3(idx),'pr')
strmax = [' ',num2str(maxy)];
text(Time3(idx),-maxy,strmax,'Color','red','FontSize',10);
ylabel('Moment, M [kips.in]]')
xlabel('Time, t [sec]')
title('Response of overturning moment of Mode 3','fontsize',F_SIZE)
 axis([0,20,-120000,120000])



subplot(4,1,4)
[maxy,idx] = max(abs(Mtotal));
plot(Time1, Mtotal, Time1(idx),Mtotal(idx),'pr')
strmax = [' ',num2str(maxy)];
text(Time1(idx),-maxy,strmax,'Color','red','FontSize',10);
ylabel('Moment, M [kips.in]')
xlabel('Time, t [sec]')
title('Response of total overturning moment','fontsize',F_SIZE)
 axis([0,20,-120000,120000])
