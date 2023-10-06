%**************************************************************************
%*     Example  9.1                                                       *
%*     filename: ch09pr01.m                                               *
%*     program listing number: 9.1                                        *
%*                                                                        *
%*     This program solves a two-dimensional non-linear equation          *
%*     by newton-raphson method.                                          *
%*                                                                        *
%*     Programed by Ryoichi Kawai for Computational Physics Course.       *
%*     Last modification:  02/04/2015.                                    *
%**************************************************************************
clear all;

% set a tolerance
tol = 1.0e-8;

% step factor (between 0 and 1)
a = 0.1;

%define functions
f1 = @(x,y) 2*x + 3*x*y - 1;
f2 = @(x,y) x*y + 3*y + 1;

%define Jacobian
J11 = @(x,y) 2+3*y;
J12 = @(x,y) 3*x;
J21 = @(x,y) y;
J22 = @(x,y) x+3;

% iniital guess
x(1) = 1;
y(1) = 0;

b=[[-f1(x(1),y(1))]; [-f2(x(1),y(1))]];
err=sqrt(b(1)*b(1)+b(2)*b(2));
if err < tol
    found = true;
else
    found = false;
end
n=1;

while not(found)
    % Construct linear equation
    A(1,1)=J11(x(n),y(n));
    A(1,2)=J12(x(n),y(n));
    A(2,1)=J21(x(n),y(n));
    A(2,2)=J22(x(n),y(n));
    % solve the linear equation
    z=A\b;
    x(n+1)=x(n)+a*z(1);
    y(n+1)=y(n)+a*z(2);
    
    % check error
    b =[[-f1(x(n+1),y(n+1))];[-f2(x(n+1),y(n+1))]];
    err=sqrt(b(1)*b(1)+b(2)*b(2));

    if err < tol
       found = true;
    else
       n=n+1;
    end
end

p=plot(x,y);
set(p,'linewidth',2,'color','black')
axis equal
xlabel(texlabel('x'),'fontsize',14)
ylabel(texlabel('y'),'fontsize',14)

xa = (-1+sqrt(7))/2;
ya = (-5+sqrt(7))/9;
fprintf('Iternations = %d\n',n)
fprintf('Solution= %d, %d\n',x(n),y(n))
fprintf('   Exact= %d, %d\n',xa,ya)
fprintf('err= %f\n',err)
