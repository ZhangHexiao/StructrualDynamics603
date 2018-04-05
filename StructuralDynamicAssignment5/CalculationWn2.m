clc;
clear;
syms k;
syms b;
M=[5/3*k-b,-2/3*k,0;-2/3*k,k-b,-1/3*k;0,-1/3*k,1/3*k-0.5*b];
f=det(M);
result=solve(f,b)

