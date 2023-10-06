%**************************************************************************
%*     Example 5.6                                                        *
%*     filename: ch05pr06.m                                               *
%*     program listing number: 5.6                                        *
%*                                                                        *
%*     This program solves Newton equation for simple harmonic oscillator *
%*     using Verlet method.                                               *
%*                                                                        *
%*     Programed by Ryoichi Kawai for Computational Physics Course        *
%*     Revised on 01/07/2014.                                             *
%**************************************************************************
clear all;

% system parameter
omega=1;

% initial conditions
x(1)=1;
v(1)=0;
t(1)=0;
x_ex(1)=cos(0);

% control parameters
tmax=8*pi/omega;
N=500;
h=tmax/N;

% the first Euler step
x(2) = x(1) + v(1)*h - omega^2*x(1)*h^2/2;

for n=2:N-1
    % Verlet method
    x(n+1)=2*x(n)-x(n-1) - omega^2*x(n)*h^2;
    t(n+1)=t(1)+n*h;
   
    % exact soution
    x_ex(n+1)=cos(omega*t(n+1));
end

% plot trajectories
subplot(1,2,1);
q=plot(t,x,'o',t,x_ex);
xlabel('t');
ylabel('displacement');
set(q(1),'Color','blue','Linewidth',2);
set(q(2),'Color','red','Linewidth',1);
legend(q,{'Verlet','Exact'});
legend(q,'Location','East');

% plot absolute error
subplot(1,2,2);
p=semilogy(t,abs(x-x_ex));
xlabel('t');
ylabel('absolute error');
set(p,'Color','blue','Linewidth',2);