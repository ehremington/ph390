%**************************************************************************
%*     Example  3.4                                                       *
%*     filename: ch03pr03.m                                               *
%*     program listing number: 3.3                                       *
%*                                                                        *
%*  This program numerically integrates x^3*exp(x)/(exp(x)-1) from        *
%*  x=0 to infinity using 8-point Gaussian Laguerre Quadrature.           *
%*                                                                        *
%*     Programed by Ryoichi Kawai for Computational Physics Course.       *
%*     Last modification:  01/06/2014.                                    *
%**************************************************************************
clear all;

N=8;

x=[1.7027963230510100e-1, 9.0370177679937991e-1,...
   2.2510866298661307,    4.2667001702876588,   ...
   7.0459054023934657,    1.0758516010180995e+1,...
   1.5740678641278005e+1, 2.2863131736889264e+1];
w=[3.6918858934163753e-1, 4.1878678081434296e-1,...
   1.7579498663717181e-1, 3.3343492261215652e-2,...
   2.7945362352256725e-3, 9.0765087733582131e-5,...
   8.4857467162725315e-7, 1.0480011748715104e-9];

for i=1:N
    f(i)=x(i)^3*exp(x(i))/(exp(x(i))-1);
end

Gauss=sum(w.*f);
Exact=pi^4/15;
fprintf('%i points Gaussian Laguerre Quadrature\n',N);
fprintf('Exact= %18.12e\n Gauss=%18.12e\n Error=%18.12e\n',...
    Exact,Gauss,abs(Exact-Gauss));
