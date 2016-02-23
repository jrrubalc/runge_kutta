function [z]=f_sys(w,t)
%
% This function calculates f_sys(w,t)
%

z=zeros(1,2);
x1=w(1);
x2=w(2);
z(1)=x2;
z(2)=(1-x1^2)*x2-x1;
%
