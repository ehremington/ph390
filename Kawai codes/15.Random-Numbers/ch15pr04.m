%**************************************************************************
%*     Example 15.4                                                       *
%*     filename: ch15pr04.m                                               *
%*     program listing number: 15.4                                       *
%*                                                                        *
%*     This program evaluates 1st through 4th moments of normal           *
%*     distribution.                                                      *
%*                                                                        *
%*     Programed by Ryoichi Kawai for Computational Physics Course.       *
%*     Last modification:  02/25/2017.                                    *
%**************************************************************************
clear all;

N=100000;

v=normrnd(0.0,1.0,[N,1]);

fprintf('order    moment\n')
for i=1:4
    m(i)=sum(v.^i)/N;
    fprintf('%3d     % 10.4e\n',i,m(i));
end
r=m(4)/m(2)^2;
fprintf('\nm4/(m2)^2 = %8.4d (exact=3)\n',r)

