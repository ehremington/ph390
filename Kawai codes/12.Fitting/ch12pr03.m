%**************************************************************************
%*     Example 12.3                                                       *
%*     filename: ch12pr03.m                                               *
%*     program listing number: 12.3                                       *
%*                                                                        *
%*     This program interpolates 11-point data with the Vandermonde       *
%*     matrix.                                                            *
%*                                                                        *
%*     Programed by Ryoichi Kawai for Computational Physics Course.       *
%*     Last modification:  02/25/2017.                                    *
%**************************************************************************
clear all

% data to be fitted
F=[0.0000, 0.6889, 0.6095, 0.0774, -0.3401, -0.3528,...
   -0.0842, 0.1620, 0.1997, 0.0681, -0.0736];
X=[0.0,1.0,2.0,3.0,4.0,5.0,6.0,7.0,8.0,9.0,10.];
N=size(F,2);

% construction of the Vandermonde matrix
x = zeros(N,N);
x(:,1)=1;
for i=2:N
    x(:,i)=X(:).^(i-1);
end

% solve the linear equation x*a=F
a=x\F';

% evaluate the function value
% between the sampling points.
M=101;
z=linspace(0,X(N),M);

for j=1:M
    y(j)=a(1);
    for i=2:N
        y(j)=y(j)+a(i)*z(j)^(i-1);
    end
end

v = sin(z).*exp(-0.2*z);

p=plot(X,F,'o',z,y,z,v,'--');
set(p,'linewidth',2);
xlabel('x','fontsize',14);
ylabel('f(x)','fontsize',14);
legend('Data','Vandermonde fit','Original');
legend('location','northeast');
