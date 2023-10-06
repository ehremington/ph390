%**************************************************************************
%*     Example  8.11                                                      *
%*     filename: ch08pr10.m                                               *
%*     program listing number: 8.10                                       *
%*                                                                        *
%*     This program solves a 2x2 linear equation by the conjugate         *
%*     gradient minimization.                                             *
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

% conjugate gradient method
n=1;
x=[1;0.5]; % starting point
u(1)=x(1);
v(1)=x(2);
r=b- A*x;
p=r;
rr=r'*r;
while abs(rr)>tol
    n=n+1;
    alpha = rr/(r'*A*r);
    x = x + alpha*p;
    u(n)=x(1);
    v(n)=x(2);
    r1=b-A*x;
    rr1=r1'*r1;
    beta=rr1/rr;
    r=r1;
    rr=rr1;
    p=r+beta*p;
end

fprintf('Solution=(%.5f,%.5f)\n',x)

p=plot(u,v);
set(p,'linewidth',2,'color','black')
axis equal tight
xlabel(texlabel('x_1'),'fontsize',14)
ylabel(texlabel('x_2'),'fontsize',14)
legend('f(x)','CG steps');
hold off
