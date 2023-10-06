%**************************************************************************
%*  Example  3.3                                                          *
%*  filename: ch03pr02.m                                                  *
%*  program listing number: 3.2                                           *
%*                                                                        *
%*  This program integrates 1/(sqrt(x)*(1+x)) from x=0 to x=1             *
%*  by removing singularity at x=0. Trapezoidal rule is used              *
%*  for the proper part of integral.                                      *
%*                                                                        *
%*     Programed by Ryoichi Kawai for Computational Physics Course.       *
%*     Last modification:  01/07/2014.                                    *
%**************************************************************************

clear all;
a = 1.0; % upper bound
N = 100; % number of segments
h = a/N; % width of segments
% integration of sqrt(x)/(1+x) with trapezoidal rule
S = sqrt(a)/(1+a)/2; %bundary value devided by 2
for i=1:N-1
x = i*h;
f = sqrt(x)/(1+x);
S = S +f;
end
proper = S*h; % integral of proper part
singular = 2*sqrt(a);
 % singular part
total = singular - proper;
exact = pi - 2*atan(1/sqrt(a));
fprintf(’Numerical = %18.12e\n’,total);
fprintf(’
 Exact = %18.12e\n’,exact);
