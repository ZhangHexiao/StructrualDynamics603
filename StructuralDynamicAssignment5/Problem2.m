
clc;
clear;
%%%%%%%Calculate the relative contribution%%%%%

Fi1 = [0.314;0.687;1];
Fi2 = [-0.5;-0.5;1];
Fi3 = [3.186;-2.186;1];
M=[1,0,0;0,1,0;0,0,0.5];

u1 = [1;2;3];
u2 = [-1;0.25;1];
u3 = [1;-1;1];

q1_1 = (Fi1')*M*u1/((Fi1')*M*Fi1);
q2_1 = (Fi2')*M*u1/((Fi2')*M*Fi2);
q3_1 = (Fi3')*M*u1/((Fi3')*M*Fi3);

q1_2 = (Fi1')*M*u2/((Fi1')*M*Fi1);
q2_2 = (Fi2')*M*u2/((Fi2')*M*Fi2);
q3_2 = (Fi3')*M*u2/((Fi3')*M*Fi3);

q1_3 = (Fi1')*M*u3/((Fi1')*M*Fi1);
q2_3 = (Fi2')*M*u3/((Fi2')*M*Fi2);
q3_3 = (Fi3')*M*u3/((Fi3')*M*Fi3);



q1_2*Fi1;
q2_2*Fi2;
q3_2*Fi3;

q1_3*Fi1
q2_3*Fi2
q3_3*Fi3