clc;
clear;

T(1) = 0;
dT = 0.1;
T1 = 2.80;
T2 = 1.28;
T3 = 0.88;

for i=1:500

    u11(i) = 0.935*cos(2.24*(T1/100*i))+0.065*cos(7.14*(T1/100*i));
    u12(i) = 2.046*cos(2.24*(T1/100*i))-0.044*cos(7.14*(T1/100*i));
    u13(i) = 2.978*cos(2.24*(T1/100*i))+0.020*cos(7.14*(T1/100*i));
    
    u21(i) = 0.105*cos(2.24*(T1/100*i))-0.438*cos(4.90*(T1/100*i))-0.667*cos(7.14*(T1/100*i));
    u22(i) = 0.230*cos(2.24*(T1/100*i))-0.438*cos(4.90*(T1/100*i))+0.458*cos(7.14*(T1/100*i));
    u23(i) = 0.334*cos(2.24*(T1/100*i))+0.875*cos(4.90*(T1/100*i))-0.209*cos(7.14*(T1/100*i));
    
    u31(i) = 0.037*cos(2.24*(T1/100*i))-0.25*cos(4.90*(T1/100*i))+1.213*cos(7.14*(T1/100*i));
    u32(i) = 0.081*cos(2.24*(T1/100*i))-0.25*cos(4.90*(T1/100*i))-0.832*cos(7.14*(T1/100*i));
    u33(i) = 0.119*cos(2.24*(T1/100*i))+0.50*cos(4.90*(T1/100*i))+0.381*cos(7.14*(T1/100*i));
    
    
    
    
    
    t(i) = i/100;
end

F_SIZE = 10;

subplot(3,1,1)

plot(t,u11)
hold on
plot(t,u12,'--b','Linewidth',1)
hold on
plot(t,u13,':r','Linewidth',2)
hold on
title('Displacement of the three-story shear building-Condition (a)','fontsize',F_SIZE+2)
ylabel('Displacement [in]')
xlabel('t/T1')

grid on;
box on

subplot(3,1,2)
plot(t,u21)
hold on
plot(t,u22,'--b','Linewidth',1)
hold on
plot(t,u23,':r','Linewidth',2)
hold on
title('Displacement of the three-story shear building-Condition (b)','fontsize',F_SIZE+2)
ylabel('Displacement [in]')
xlabel('t/T1')


grid on;
box on

subplot(3,1,3)
plot(t,u31)
hold on
plot(t,u32,'--b','Linewidth',1)
hold on
plot(t,u33,':r','Linewidth',2)
hold on
title('Displacement of the three-story shear building-Condition (c)','fontsize',F_SIZE+2)
ylabel('Displacement [in]')
xlabel('t/T1')

h_legend=legend('U1','U2','U3');
grid on;
box on








