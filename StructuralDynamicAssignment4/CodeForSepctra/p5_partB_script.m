clear all;
clc;
LP=load('LomaPrieta.th');

p=386.1*transpose(LP(1:6001));
t=0:0.005:30;
Tn1=0.5;
m=200/386.1;
[u1,v,acc,ab_acc,fs1,k1]=centraldifferencesolver(Tn1,0.005,p,0.05,m);
ry1=2;
u01=max(abs(u1))
f01=max(abs(fs1))
fy1=f01/ry1
uy1=fy1/k1
[un1,vn,accn,ab_accn,fsn1]=nonlinearsolver(Tn1,0.005,p,0.05,fy1,m);
um1=max(abs(un1))
miu1=um1/uy1
x1=un1(1:end-1);
subplot(4,2,1);
plot(t,un1)
ylabel('response u [in]')
title('Tn=0.5; Ry=2')
subplot(4,2,2);
plot(x1,fsn1)
ylabel('resisting force [kips]')
title('Tn=0.5; Ry=2')

Tn2=0.5;
m=200/386.1;
[u2,v,acc,ab_acc,fs2,k2]=centraldifferencesolver(Tn2,0.005,p,0.05,m);
ry2=4;
u02=max(abs(u2))
f02=max(abs(fs2))
fy2=f02/ry2
uy2=fy2/k2
[un2,vn,accn,ab_accn,fsn2]=nonlinearsolver(Tn2,0.005,p,0.05,fy2,m);
um2=max(abs(un2))
miu2=um2/uy2
x2=un2(1:end-1);
subplot(4,2,3);
plot(t,un2)
ylabel('response u [in]')
title('Tn=0.5; Ry=4')
subplot(4,2,4);
plot(x2,fsn2)
ylabel('resisting force [kips]')
title('Tn=0.5; Ry=4')

Tn3=1;
m=200/386.1;
[u3,v,acc,ab_acc,fs3,k3]=centraldifferencesolver(Tn3,0.005,p,0.05,m);
ry3=2;
u03=max(abs(u3))
f03=max(abs(fs3))
fy3=f03/ry3
uy3=fy3/k3
[un3,vn,accn,ab_accn,fsn3]=nonlinearsolver(Tn3,0.005,p,0.05,fy3,m);
um3=max(abs(un3))
miu3=um3/uy3
x3=un3(1:end-1);
subplot(4,2,5);
plot(t,un3)
ylabel('response u [in]')
title('Tn=1.0; Ry=2')
subplot(4,2,6);
plot(x3,fsn3)
ylabel('resisting force [kips]')
title('Tn=1.0; Ry=2')

Tn4=1;
m=200/386.1;
[u4,v,acc,ab_acc,fs4,k4]=centraldifferencesolver(Tn4,0.005,p,0.05,m);
ry4=4;
u04=max(abs(u4))
f04=max(abs(fs4))
fy4=f04/ry4
uy4=fy4/k4
[un4,vn,accn,ab_accn,fsn4]=nonlinearsolver(Tn4,0.005,p,0.05,fy4,m);
um4=max(abs(un4))
miu4=um4/uy4
x4=un4(1:end-1);
subplot(4,2,7);
plot(t,un4)
xlabel('Time t [s]')
ylabel('response u [in]')
title('Tn=1.0; Ry=4')
subplot(4,2,8);
plot(x4,fsn4)
xlabel('Displacement u [in]')
ylabel('resisting force [kips]')
title('Tn=1.0; Ry=4')
