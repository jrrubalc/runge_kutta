%
% This code implements Dormand & Prince's fifth order method (a six
% stage fifth order Runge Kutta method) to solve an ODE system. 
% After the calculation, it saves the workspace to a data file.
%
clear
%
m=2;
w0=[1, 0];
h=0.1;
nstep=40/h;
%
w=zeros(nstep+1,m);
t=zeros(nstep+1,1);
t(1)=0;
w(1,1:m)=w0;
%
p=6;
d=[0, 1/5, 3/10, 4/5, 8/9, 1];
c=[0, 0, 0, 0, 0, 0;
   1/5, 0, 0, 0, 0, 0;
   3/40, 9/40, 0, 0, 0, 0;
   44/45, -56/15, 32/9, 0, 0, 0;
   19372/6561, -25360/2187, 64448/6561, -212/729, 0, 0;
   9017/3168, -355/33, 46732/5247, 49/176, -5103/18656, 0];
b=[35/384, 0, 500/1113, 125/192, -2187/6784, 11/84];
k=zeros(p,m);
%
for j=1:nstep,
  for i=1:p,
    k(i,1:m)=h*f_sys(w(j,1:m)+c(i,1:i-1)*k(1:i-1,1:m), t(j)+d(i)*h);
  end
  w(j+1,1:m)=w(j,1:m)+b*k;
  t(j+1)=t(j)+h;
end
%
save data_sRK5
%
