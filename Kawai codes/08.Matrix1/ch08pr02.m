%**************************************************************************
%*     Example  8.2                                                       *
%*     filename: ch08pr02.m                                               *
%*     program listing number: 8.2                                        *
%*                                                                        *
%*     This program solves a upper-triangular linear equation with        *
%*     the backsubstitution method.                                       *
%*                                                                        *
%*     Programed by Ryoichi Kawai for Computational Physics Course.       *
%*     Last modification:  10/13/2013.                                    *
%**************************************************************************
clear all;

% define matrix A and vector b
A=[[3, -1, 4];[0,2,-1];[0,0,2]];
b=[-1;-2;4];

% backsubstitution
for i=3:-1:1
    Ax=0;
    for j=i+1:3
        Ax = Ax+A(i,j)*x(j);
    end
    x(i) = (b(i)-Ax)/A(i,i);
end

fprintf('x=%.1f, y=%.1f, z=%.1f\n',x)

