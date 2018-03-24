close all
clear all
clc

% 
% x=0:0.01:5;
% 
% ForcedVibration = 1 - (x.*2.*pi).^(-1).*sin(x.*2.*pi);
% 
% y=0:0.01:5;
% a=(2.*y.*pi).^(-1);
% c=(1-a.*(sin(2.*pi.*y))).^2;  
% b=(pi*pi*y.*y).^(-1).*(sin(pi*y)).^4;
% FreeVibration = (c+b).^0.5;
% 
% 
% 
% grid on
% hold on
% plot(y,FreeVibration,'--');
% 
% plot(x,ForcedVibration);
% xlabel('to/Tn');
% ylabel('Rd')



y=0:0.01:5;
a=(2.*y.*pi).^(-1);
c=(1-a.*(sin(2.*pi.*y))).^2;  
b=(pi*pi*y.*y).^(-1).*(sin(pi*y)).^4;
FreeVibration = (c+b).^0.5;


plot(y,FreeVibration);
xlabel('to/Tn');
ylabel('Rd')
grid on






% x1=0:0.01:0.5;
% y1 = 2.*x1-1/pi.*sin(2.*pi.*x1);
% 
% x2=0.5:0.01:5;
% y2 = cos(2*pi*x2-pi)+1/pi*sin(2*pi.*x2-pi)-1/pi*sin(2*pi.*x2) ;
% plot(x1,y1,x2,y2);
% grid on
% xlabel('t/Tn')
% ylabel('u(t)/(ust)0')


% x3=0:0.01:2;
% y3=0.5*x3-0.25/pi.*sin(2.*pi.*x3);
% 
% x4=2:0.01:5;
% y4 = cos(2*pi*x4-4*pi)+0.25/pi*sin(2*pi.*x4-4*pi)-0.25/pi*sin(2*pi.*x4) ;
% plot(x3,y3,x4,y4);
% xlabel('t/Tn')
% ylabel('u(t)/(ust)0')
% grid on
% 











