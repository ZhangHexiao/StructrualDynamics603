clc;
clear;
%%%%%Orthogonality%%%%%
Fi1 = [0.314;0.687;1];
Fi2 = [-0.5;-0.5;1];
Fi3 = [3.186;-2.186;1];
M=[1,0,0;0,1,0;0,0,0.5];
Fi2'*M*Fi3;
%%%%%%%Normalized Modes%%%%%%%
M1=Fi1'*M*Fi1;
M2=Fi2'*M*Fi2;
M3=Fi3'*M*Fi3
Fi1/(M1^0.5);
Fi2/(M2^0.5);
Fi3/(M3^0.5)
