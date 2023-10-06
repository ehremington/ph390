%**************************************************************************
%*     Example 12.5                                                       *
%*     filename: ch12pr05.m                                               *
%*     program listing number: 12.5                                       *
%*                                                                        *
%*     This program interpolates 11-point data with the linear            *
%*     regression method.                                                 *
%*                                                                        *
%*     Programed by Ryoichi Kawai for Computational Physics Course.       *
%*     Last modification:  02/25/2017.                                    *
%**************************************************************************
clear all;

% Data set (no error bar)
x=[0.0,1.0,2.0,3.0,4.0,5.0,6.0,7.0,8.0,9.0,10.];
y=[0.1,0.90,1.7,3.4,4.5,4.7,6.2,7.6,7.85,9.03,9.6];
N=size(y,2);

% Linear regression
F=sum(y);
X=sum(x);
X2=sum(x.^2);
XF=sum(x.*y);
b=(F*X2-X*XF)/(N*X2-X^2);
a=(N*XF-X*F)/(N*X2-X^2);

% fitted curve
f=a*x+b;

p=plot(x,f);
set(p,'linewidth',2,'color','red')
hold on
r=plot(x,y,'o');
set(r,'linewidth',2,'color','black')
xlabel('x','fontsize',14)
ylabel('f(x)','fontsize',14)
hold off
