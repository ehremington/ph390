%**************************************************************************
%*     Example 5.4                                                        *
%*     filename: ch05pr04.m                                               *
%*     program listing number: 5.4                                        *
%*                                                                        *
%*     This program solves Newton equation for interacting two cars       *
%*     using Runge-Kutta 2nd order methods.                               *
%*                                                                        *
%*     Programed by Ryoichi Kawai for Computational Physics Course        *
%*     Revised on 01/07/2014.                                             *
%**************************************************************************
clear all;

% Control parameters
tmax=8; N=400; h=tmax/N;

% initial conditions
v1(1)=1.2; v2(1)=1.0; t(1)=0;

% 2nd-order Runge-Kutta method
for n=1:N-1
    t(n+1) = t(1)+n*h;
    k1 =   v2(n)-v1(n);
    l1 = -(v2(n)-v1(n));
    mid1 = v1(n)+k1*h/2;
    mid2 = v2(n)+l1*h/2;
    k2 =   mid2-mid1;
    l2 = -(mid2-mid1);
    v1(n+1)=v1(n)+k2*h;
    v2(n+1)=v2(n)+l2*h;
end

subplot(1,2,1); 
p=plot(t,v1,t,v2);
xlabel('t');
ylabel(texlabel('velocity'));
set(p(1),'Color','blue','Linewidth',2);
set(p(2),'Color','red','Linewidth',2);
legend(p,{texlabel('v_1'),texlabel('v_2')});
legend(p,'Location','SouthEast');

subplot(1,2,2);
q=plot(t,v1-v2);
xlabel('t');
ylabel(texlabel('v_1-v_2'));
set(q,'Linewidth',2);
