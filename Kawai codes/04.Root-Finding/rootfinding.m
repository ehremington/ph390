%**************************************************************************
%*     Secion 4.6.                                                        *
%*     filename: rootfinding.m                                            *
%*     program listing number: 4.6-2                                      *
%*                                                                        *
%*     Inputs                                                             *
%*         f = function name                                              *
%*         x1, x2 = bracket                                               *
%*         N = interations of bisection                                   *
%*         tol = tolerance for scant method                               *
%*     Output:                                                            *
%*         x = root                                                       *
%*                                                                        *
%*     This program find a root of a given funcion f(x)=0                 *
%*     using the bisection and secant root finding method                 *
%*     finding methods.                                                   *
%*                                                                        *
%*     Programed by Ryoichi Kawai for Computational Physics Course        *
%*     Revised on 01/07/2014.                                             *
%**************************************************************************

function x=rootfinding(f,x1,x2,N,tol)

% check if the initial bracket satisfies the necessary condition.
f1=f(x1);
f2=f(x2);
if f1*f2 > 0
   error('Bracket is incorrect');
   display(x1,x2,f1,f2);
end

% Bisection method (N iteration)
n=0; % iteration counter

% mid point
xm = (x1+x2)/2;
fm = f(xm);
n=0;

while n<N
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

% Secant method
dx = (x2-x1)/10;
x1 = xm;
f1 = fm;
x2 = x1 + dx;
f2 = f(x2);
n = 0;
while abs(f2)> tol
    x = x2 - (x2-x1)/(f2-f1)*f2;
    x1 = x2;
    f1 = f2;
    x2 = x;
    f2 = f(x);
    n = n+1;
end
