%**************************************************************************
%*     Example 12.2                                                       *
%*     filename: ch12pr02.m                                               *
%*     program listing number: 12.2                                       *
%*                                                                        *
%*     This program interpolates 11-point data with cubic spline.         *
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

% Set up linear equation A*P=G
for j=1:N-1
    G(j)=3*((F(j+2)-F(j+1))/h(j+1)-(F(j+1)-F(j))/h(j));
end

A=zeros(N-2,N-2);
for j=1:N-1
    A(j,j)=(h(j+1)+h(j))/2;
end
for j=1:N-2
    A(j,j+1)=h(j+1);
    A(j+1,j)=h(j+1);
end

P=A\G'; %Solve the linear equation

% Shift the index to our convention.
for j=1:N-1
    Q(j+1)=P(j);
end

% Amend the boundary values.
Q(1)=0;
Q(N+1)=0;

M=10;
dt=1.0/M;
t=linspace(0.,dt*(M-1),M);

n=0;
for k=1:N
    for j=1:M
        n=n+1;
        x(n)=X(k)+t(j)*h(k);
        y(n)=h(k)^2/6 * Q(k)   * t(j)*(t(j)+1)*(t(j)-1) ...
            -h(k)^2/6 * Q(k+1) * t(j)*(t(j)-1)*(t(j)-2) ...
            +F(k+1)*t(j) + (1-t(j))*F(k);
    end
end

v=sin(x).*exp(-0.2*x);
p=plot(x,y,X,F,'o',x,v,'--');
set(p,'linewidth',2);
xlabel('x','fontsize',14);
ylabel('f(x)','fontsize',14);
legend('Cubic Spline','Data ', 'Original');
legend('location','northeast');
