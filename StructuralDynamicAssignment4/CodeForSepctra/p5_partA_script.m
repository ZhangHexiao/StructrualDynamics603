clear all;
clc;
LP=load('LomaPrieta.th');

p=386.1*transpose(LP(1:6001));
Tn1=0.5;
m=200/386.1;
t=0:0.005:30;
[u1,v,acc,ab_acc,fs,k1]=centraldifferencesolver(Tn1,0.005,p,0.05,m);
u01=max(abs(u1))
f01=u01*k1
plot(t,u1)
hold on 
Tn2=1;
[u2,v,acc,ab_acc,fs,k2]=centraldifferencesolver(Tn2,0.005,p,0.05,m);
u02=max(abs(u2))
f02=u02*k2
plot(t,u2)
xlabel('Time t [s]')
ylabel('Linear Response u [in]')
legend('Tn=0.5','Tn=1.0')
hold off
