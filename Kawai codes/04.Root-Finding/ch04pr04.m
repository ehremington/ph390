%**************************************************************************
%*     Example 4.4                                                        *
%*     filename: ch04pr04.m                                               *
%*     program listing number: 4.4                                        *
%*                                                                        *
%*     This program seeks the root of                                     *
%*             cos(3*x)*sin(x)=0                                          *
%*     between x=0.2 and 0.8 using bisection, Newton-Raphson, and         *
%*     secant methods.                                                    *
%*                                                                        *
%*     Programed by Ryoichi Kawai for Computational Physics Course        *
%*     Revised on 01/07/2014.                                             *
%**************************************************************************
clear all;

% define function and its derivative
f=@(x) cos(3*x)*sin(x);
df=@(x) -3*sin(3*x)*sin(x) + cos(3*x)*cos(x);

% initial bracket
x1=0.2;
x2=0.8;
f1=f(x1);
f2=f(x2);
if f1*f2 > 0
   error('Bracket is incorrect');
end

% tolerance
epsilon=1e-8;

% bisection 10 iterations

% iteration counter
n=0;

% Bisection method (10 iterations)
xm = (x1+x2)/2;
fm = f(xm);
while n<10
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

fprintf('Bisection      = %.8f  (iteration= %i)\n',xm,n);

% Newton-Raphson method
x = xm;
fx = fm;
n = 0;
while abs(fx)> epsilon
    dfx = df(x);
    x = x - fx/dfx;
    fx = f(x);
    n = n+1;
end

fprintf('Newton-Raphson = %.8f  (iteration= %i)\n',x,n);

% Secant method
dx = (x2-x1)/10;
x1 = xm;
f1 = fm;
x2 = x1 + dx;
f2 = f(x2);
n = 0;
while abs(f2)> epsilon
    x = x2 - (x2-x1)/(f2-f1)*f2;
    x1 = x2;
    f1 = f2;
    x2 = x;
    f2 = f(x);
    n = n+1;
end

fprintf('Secant         = %.8f  (iteration= %i)\n',x,n);

