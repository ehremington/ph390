%**************************************************************************
%*     Example  8.10                                                      *
%*     filename: ch08pr09.m                                               *
%*     program listing number: 8.9                                        *
%*                                                                        *
%*     This program solves a 2x2 linear equation by the steepest descent  *
%*     minimization.                                                      *
%*                                                                        *
%*     Programed by Ryoichi Kawai for Computational Physics Course.       *
%*     Last modification:  10/13/2013.                                    *
%**************************************************************************
clear all;
A=[[4,1];[1,3]];
b=[1;2];
c=[-0.65,  -0.5, -0.3,  -0.1, 0.1,  0.3, 0.5, 0.7, 0.9, 1.1];
tol = 1e-8;

% contour plot of the cost function
x=linspace(-1,1.2);
y=linspace(-1,2);
[X,Y]=meshgrid(x,y);
N=size(X,2);
M=size(Y,2);
for i=1:N
    for j=1:M
        Z(i,j) = 0.5*(X(i,j)^2*A(1,1)+(A(1,2)+A(2,1))*X(i,j)*Y(i,j)...
                 +A(2,2)*Y(i,j)^2) - X(i,j)*b(1)-Y(i,j)*b(2);
    end
end

contour(X,Y,Z,c);
hold on

% steepest descent with line minimization
n=1;
x=[1;0.5]; % starting point
u(1)=x(1);
v(1)=x(2);
g = A*x-b;
gg=g'*g;
while abs(gg)>tol
    n=n+1;
    lambda = gg/(g'*A*g);
    x = x - lambda * g;
    u(n)=x(1);
    v(n)=x(2);
    g = A*x-b;
    gg = g'*g;
end

fprintf('Solution=(%.5f,%.5f)\n',x)

p=plot(u,v);
set(p,'linewidth',2,'color','black')
axis equal tight
xlabel(texlabel('x_1'),'fontsize',14)
ylabel(texlabel('x_2'),'fontsize',14)

hold off
