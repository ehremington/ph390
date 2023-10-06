%**************************************************************************
%*     Example  8.1                                                       *
%*     filename: ch08pr01.m                                               *
%*     program listing number: 8.1                                        *
%*                                                                        *
%*     This program checks the properties of triangular matrices.         *
%*                                                                        *
%*     Programed by Ryoichi Kawai for Computational Physics Course.       *
%*     Last modification:  01/31/2015.                                    *
%**************************************************************************
clear all;

% define matrices A and B
A=[[2, 0, 0];[-1,1,0];[3,2,-1]]
B=[[1, 0, 0];[2,4,0];[-1,-2,3]]

C=A*B;
fprintf('A*B\n')
fprintf('%3d %3d %3d\n',C')
fprintf('\nInverse of A\n')
D=inv(A);
fprintf('%15.4e %15.4e %15.4e\n',D')
E2=det(A);
fprintf('Determinant of A = %d\n',E2)


