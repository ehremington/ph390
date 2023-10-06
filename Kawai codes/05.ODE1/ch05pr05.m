%**************************************************************************
%*     Example 5.5                                                        *
%*     filename: ch05pr05.m                                               *
%*     program listing number: 5.5                                        *
%*                                                                        *
%*     This program solves Newton equation for simple harmonic oscillator *
%*     using Runge-Kutta 4th order methods.                               *
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
h=tmax/N

for n=1:N-1
    % 4th-order Runge-Kutta
    kv1=-omega^2*x(n);
    kx1=v(n);
    
    v_mid = v(n)+kv1*h/2;
    x_mid = x(n)+kx1*h/2;
    kv2 = -omega^2*x_mid;
    kx2 = v_mid;
      
    v_mid = v(n)+kv2*h/2;
    x_mid = x(n)+kx2*h/2;
    kv3 = -omega^2*x_mid;
    kx3 = v_mid;

    v_end = v(n)+kv3*h;
    x_end = x(n)+kx3*h;
    kv4 = -omega^2*x_end;
    kx4 = v_end;
    
    v(n+1)=v(n)+(kv1+2*(kv2+kv3)+kv4)*h/6;
    x(n+1)=x(n)+(kx1+2*(kx2+kx3)+kx4)*h/6;
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
legend(q,{'RK4','Exact'});
legend(q,'Location','East');

% plot absolute error
subplot(1,2,2);
p=semilogy(t,abs(x-x_ex));
xlabel('t');
ylabel('absolute error');
set(p,'Color','blue','Linewidth',2);

