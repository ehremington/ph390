%**************************************************************************
%*     Section 5.5.1                                                      *
%*     filename: ch05pr07.m                                               *
%*     program listing number: 5.7                                        *
%*                                                                        *
%*     This program solves the Brusselator model                          *
%*     using Runge-Kutta 4th order methods.                               *
%*                                                                        *
%*     Programed by Ryoichi Kawai for Computational Physics Course        *
%*     Revised on 01/07/2014.                                             *
%**************************************************************************
clear all;

% fixed parameter
a=1;

% control parameter
b=input('Enter value for b [1.5-2.5] =');

%initial conditions
x(1)=1;
y(1)=1;
t(1)=0;

% duration 
tmax=100;

% number of integration steps
N=2000;

% step size
h=tmax/N;

for n=1:N-1
    % 4th-order Runge-Kutta
    kx1=a-(b+1)*x(n)+x(n)^2*y(n);
    ky1=b*x(n)-x(n)^2*y(n);
    
    x_mid = x(n)+kx1*h/2;
    y_mid = y(n)+ky1*h/2;
    kx2 = a - (b+1)*x_mid+x_mid^2*y_mid;
    ky2 = b*x_mid-x_mid^2*y_mid;
      
    x_mid = x(n)+kx2*h/2;
    y_mid = y(n)+ky2*h/2;
    kx3 = a - (b+1)*x_mid+x_mid^2*y_mid;
    ky3 = b*x_mid-x_mid^2*y_mid;

    x_end = x(n)+kx3*h;
    y_end = y(n)+ky3*h;
    kx4 = a - (b+1)*x_end+x_end^2*y_end;
    ky4 = b*x_end-x_end^2*y_end;
    
    x(n+1)=x(n)+(kx1+2*(kx2+kx3)+kx4)*h/6;
    y(n+1)=y(n)+(ky1+2*(ky2+ky3)+ky4)*h/6;
    t(n+1)=t(1)+n*h;
end

% plot individual trajectories
subplot(1,2,1);
p=plot(t,x,t,y);
xlabel('t');
ylabel('concentration');
set(p(1),'Color','blue','Linewidth',2);
set(p(2),'Color','red','Linewidth',2);
legend(p,{'x','y'});
legend(p,'Location','SouthEast');

% plot phase trajectory
subplot(1,2,2);
q=plot(x,y);
xlabel('x');
ylabel('y');
set(q(1),'color','blue');