%**************************************************************************
%*     Example  8.6                                                       *
%*     filename: ch08pr06.m                                               *
%*     program listing number: 8.6                                        *
%*                                                                        *
%*     This program solves a simple linear equation with LU decomposition.*
%*     MATLAB function lu() is used.                                      *
%*                                                                        *
%*     Programed by Ryoichi Kawai for Computational Physics Course.       *
%*     Last modification:  10/13/2013.                                    *
%**************************************************************************
clear all;

% Define a matrix
A=[[3, -1, 4];[2, 0, -1];[0, 3, 2]];
b=[2;-1;3];

% LU dcomposition
[L U P]=lu(A);

% Rcover the original matrix
S=P*L*U;

% Show the results
fprintf('\nA (Original Matrix)\n')
fprintf('%8.5f  %8.5f  %8.5f\n',A')
fprintf('\nL (Lower Triangular Matrix)\n')
fprintf('%8.5f  %8.5f  %8.5f\n',L')
fprintf('\nU (Upper Triangular Matrix)\n')
fprintf('%8.5f  %8.5f  %8.5f\n',U')
fprintf('\nP (Permutation Matrix)\n')
fprintf('%i  %i  %i\n',P')
fprintf('\nP*L*U\n')
fprintf('%8.5f  %8.5f  %8.5f\n',S')

b = P*b;
% forward substition
for i=1:3
    Ly=0;
    for j=1:i-1
        Ly = Ly+L(i,j)*y(j);
    end
    y(i) = (b(i)-Ly)/L(i,i);
end

% backsubstitution
for i=3:-1:1
    Ux=0;
    for j=i+1:3
        Ux = Ux+U(i,j)*x(j);
    end
    x(i) = (y(i)-Ux)/U(i,i);
end

fprintf('\nx=%.5f,  y=%.5f, z=%.5f\n',x)




