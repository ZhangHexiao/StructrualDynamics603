clc;
clear;

m = 100/386*[1,0,0;0,1,0;0,0,0.5];
k = 168/9*[16,-7,0;-7,10,-3;0,-3,3];

a = [1/12.01,12.01;1/38.90,38.90];
cita = [0.05;0.05];
b=2*a^(-1)*cita;
b(1,1);
b(2,1);
c=b(1,1)*m + b(2,1)*k;

Fi = [0.6375,0.9827,1.5778;1.2750,0.9829,-1.1270;1.9125,-1.9642,0.4508];

M = Fi'*m*Fi;

K = Fi'*k*Fi;

C = Fi'*c*Fi
