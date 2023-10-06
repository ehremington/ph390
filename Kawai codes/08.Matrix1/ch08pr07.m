%**************************************************************************
%*     Example  8.7                                                       *
%*     filename: ch08pr07.m                                               *
%*     program listing number: 8.7                                        *
%*                                                                        *
%*     This program solves a tridiagonal system with backward elimination.*
%*     Then, find solution by forward substitution.                       *
%*                                                                        *
%*     Programed by Ryoichi Kawai for Computational Physics Course.       *
%*     Last modification:  02/01/2015.                                    *
%**************************************************************************
clear all

% Define matrices. No need to use the full matrix.
d=[2,3,4,3]; % diagonal elements
u=[2,3,3,0]; % above diagonal
l=[0,2,3,3]; % below diagonal
b=[1,2,3,4]; % right hand side

% Calculation of determinant
D(1)=d(1);
D(2)=d(2)*D(1)-l(1)*u(1);
for i=3:4
    D(i)=d(i)*D(i-1)-l(i-1)*u(i-1)*D(i-2);
end

fprintf('Determinant %d\n',D(4))
if D(4) == 0
    fprintf('Singular')
    stop
end

% Decomposition by backword elimination
Y(3)=-l(4)/d(4);
Z(3)= b(4)/d(4);
for i=3:-1:2
    Y(i-1)=-l(i)/(d(i)+u(i)*Y(i));
    Z(i-1)=(b(i)-u(i)*Z(i))/(d(i)+u(i)*Y(i));
end

% Forward substitution
x(1)=(b(1)-u(1)*Z(1))/(d(1)+u(1)*Y(1));
for i=1:3
    x(i+1)=Y(i)*x(i)+Z(i);
end

% Answer
fprintf('x= %f %f %f %f\n',x)

% Check the errors. 
s(1)=d(1)*x(1)+u(1)*x(2)-b(1);
s(2)=l(2)*x(1)+d(2)*x(2)+u(2)*x(3)-b(2);
s(3)=l(3)*x(2)+d(3)*x(3)+u(3)*x(4)-b(3);
s(4)=l(4)*x(3)+d(4)*x(4)-b(4);
fprintf('Error= %e %e %e %e\n',s)
