%**************************************************************************
%*     Example 12.6                                                       *
%*     filename: ch12pr06.m                                               *
%*     program listing number: 12.6                                       *
%*                                                                        *
%*     This program interpolates 11-point data with the quadratic         *
%*     regression method.                                                 *
%*                                                                        *
%*     Programed by Ryoichi Kawai for Computational Physics Course.       *
%*     Last modification:  02/25/2017.                                    *
%**************************************************************************
clear all

% Generate a noisy data set
N=13;
sm=5.0;
for i=1:N
    x(i)=i-7.+random('unif',-0.2,0.2);
    s(i)=random('norm',0,sm/2)+sm;
    y(i)=-2*x(i)^2+s(i)*random('unif',0.2,0.9);
end

M=3; % number of parameters

% Construct Jacobian matrix and a vector
for i=1:N
    for j=1:M
        J(i,j)=x(i)^(j-1)/s(i);
    end
    b(i)=y(i)/s(i);
end

% Solve the linear equation
A=J'*J;
c=J'*b';
lambda=A\c;

% constructe the fitted curve
K=121;
z=linspace(-6.0,6.0,K);
f=lambda(1)+lambda(2)*z+lambda(3)*z.^2;

r=errorbar(x,y,s,'o');
set(r,'linewidth',2,'color','black')
hold on

p=plot(z,f);
set(p,'linewidth',2,'color','red')
xlabel('x','fontsize',14)
ylabel('f(x)','fontsize',14)
legend('data','fitted curve')
legend('location','south')
xlim([-7 7]);
hold off

