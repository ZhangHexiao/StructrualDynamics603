clc;
clear;
syms k;
syms b;

% b=0.209;
%    b =1;
 b=2.124;

M=[5/3-b,-2/3,0;-2/3,1-b,-1/3;0,-1/3,1/3-b*0.5]
syms fi1;
syms fi2;
Fi1 = [fi1;fi2;1];
x=M*Fi1;
x1 = x(1,1)
x2 = x(2,1)
x3 = x(3,1)
solve(x1,x2,fi1,fi2)

