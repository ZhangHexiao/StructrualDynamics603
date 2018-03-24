%Problem 1 - Central Difference Method Solver
function [u,v,acc,ab_acc,fs,k]=centraldifferencesolver(Tn,dt,accg,z,m)


wn=2*pi/Tn;
k=m*wn^2;
c=z*2*m*wn;
u0=0;u(1)=u0;
v0=0;v(1)=v0;
p=-m*accg;
acc(1)=(p(1)-c*v0-k*u0)/m;
kh=m/dt^2+c/(2*dt);
alpha=m/dt^2-c/(2*dt);
beta=k-2*m/dt^2;
u_1 = u0-dt*v0+dt^2/2*acc(1);
ph(1) = p(1)-alpha*u_1-beta*u(1);

for i = 2:(length(p)-1)    

u(i) = ph(i-1)/kh;
ph(i) = p(i)-alpha*u(i-1)-beta*u(i);
u(i+1) = ph(i)/kh;
v(i) = (u(i+1)-u(i-1))/(2*dt);
acc(i) = (u(i+1)-2*u(i)+u(i-1))/(dt^2);
ab_acc(i) =acc(i)+accg(i);
end
fs=k*u;
fd=c*v;

end



    