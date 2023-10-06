%**************************************************************************
%*     Example  4.3                                                       *
%*     filename: ch04pr03.m                                               *
%*     program listing number: 4.3                                        *
%*                                                                        *
%*     This program finds roots of a cubic equation                       *
%*               a*x^3 + b*x^2 + c*x + d = 0                              *
%*     using the bisection method.                                        *
%*                                                                        *
%*     Programed by Ryoichi Kawai for Computational Physics Course.       *
%*     Last modification:  01/24/2018.                                    *
%**************************************************************************
clear all;

% define a cubic equation
b = -9;
c = 23;
d =-15;

f = @(x) x^3+b*x^2+c*x+d;

% tolerance
epsilon=input('Enter tolerance =');

% initial bracket
x1=1.5;
x2=4;
f1 = f(x1);
f2 = f(x2);

% iteration counter
n=0;

% mid point
xm = (x1+x2)/2;
fm = f(xm);
while x2-x1 > epsilon
    if f1*fm < 0  % root in the lower half
        x2=xm;
        f2=fm;
    else          % root in the upper half
        x1=xm;
        f1=fm;
    end
    xm = (x1+x2)/2;  % new mid point
    fm = f(xm);
    n=n+1;
end

fprintf('Answer = %.5f  (iteration= %i)\n',xm,n);
