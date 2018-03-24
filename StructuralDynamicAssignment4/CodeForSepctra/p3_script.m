clear all;
clc;
Canoga=load('Canoga_Park.th');
Canoga2=[Canoga;0; 0; 0; 0; 0; 0];
p=386.1*transpose(Canoga2);

for f=1:40
    Tn=f*0.2;
    
[ul,vl,accl,ab_accl,fsl,kl]=centraldifferencesolver(Tn,0.01,p,0.05,1);
u0=max(abs(ul));
f0=max(abs(fsl));
for i=1:1000
    fybar(i)=0.001*i;
    fy(i)=fybar(i)*f0;
    uy(i)=fy(i)/kl;
[un,vn,accn,ab_accn,fsn]=nonlinearsolver(Tn,0.01,p,0.05,fy(i),1);
um(i)=max(abs(un));
miu(i) =um(i)/uy(i);
end
xq=0.1:0.1:1000;
miuinterp=interp1(miu,xq);
fybarinterp=interp1(fybar,xq);
k(1)=find(abs(miuinterp-1)<0.01,1,'first');
k(2)=find(abs(miuinterp-1.5)<0.01,1,'first');
k(3)=find(abs(miuinterp-2)<0.01,1,'first');
k(4)=find(abs(miuinterp-4)<0.01,1,'first');
k(5)=find(abs(miuinterp-8)<0.01,1,'first');
for n=1:5
    fy_bar(n)=fybarinterp(k(n));
end 
for j=1:5
    f_y(j)=f0*fy_bar(j);
[un,vn,accn,ab_accn,fsn]=nonlinearsolver(Tn,0.01,p,0.05,f_y(j),1);
ab_accmax(j)=max(abs(ab_accn));
end 
amiu1(f)=ab_accmax(1);
amiu15(f)=ab_accmax(2);
amiu2(f)=ab_accmax(3);
amiu4(f)=ab_accmax(4);
amiu8(f)=ab_accmax(5);
end 


figure
F_SIZE = 10;
tn=0.2:0.2:8;
plot(tn,amiu1,'Linewidth',2)
hold on
plot(tn,amiu15,'Linewidth',2.5)
hold on
plot(tn,amiu2,'Linewidth',3)
hold on
plot(tn,amiu4,'Linewidth',3.5)
hold on
plot(tn,amiu8,'Linewidth',4)
hold on
legend('u=1','u=1.5','u=2','u=4','u=8');
ylabel('Absolute Acceleration [in/sec^{2}]')
xlabel('Period, t [sec]')
title('Constant ductility absolute acceleration response spectra','fontsize',F_SIZE+2)
grid on;
box on












