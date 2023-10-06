%**************************************************************************
%*     Example 12.1                                                       *
%*     filename: ch12pr01.m                                               *
%*     program listing number: 12.1                                       *
%*                                                                        *
%*     This program interpolates 11-point data with linear spline.        *
%*                                                                        *
%*     Programed by Ryoichi Kawai for Computational Physics Course.       *
%*     Last modification:  02/25/2017.                                    *
%**************************************************************************
clear all

% data to be fitted
F=[0.0000, 0.6889, 0.6095, 0.0774, -0.3401, -0.3528,...
   -0.0842, 0.1620, 0.1997, 0.0681, -0.0736];
X=[0.0,1.0,2.0,3.0,4.0,5.0,6.0,7.0,8.0,9.0,10.];
N=size(F,2)-1;
h=X(2:N+1)-X(1:N);

M=10;
dt=1.0/M;
t=linspace(0.,dt*(M-1),M);

n=0;
for i=1:N-1
    % linear interpolation between two adjacent data points
    for j=1:M
        n=n+1;
        x(n)=t(j)*h(i)+X(i);
        y(n)=(1-t(j))*F(i)+t(j)*F(i+1);
    end
end

n=n+1;
x(n)=X(N+1);
y(n)=F(N+1);
z=sin(x).*exp(-0.2*x);

p=plot(x,y,X,F,'o',x,z,'--');
set(p,'linewidth',2);
xlabel('x','fontsize',14);
ylabel('f(x)','fontsize',14);
legend('Spline','Data','Original');
legend('location','northeast');
