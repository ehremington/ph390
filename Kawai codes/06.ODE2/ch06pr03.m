%**************************************************************************
%*     Section  6.3.1                                                     *
%*     filename: ch06pr03.m                                               *
%*     program listing number: 6.3                                        *
%*                                                                        *
%*     This program finds the wave function of freely falling particle.   *
%*     to reach height yf in travel time tf.                              *
%*                                                                        *
%*     Programed by Ryoichi Kawai for Computational Physics Course.       *
%*     Last modification:  10/13/2013.                                    *
%**************************************************************************
clear all;
clc

% define w(x) in Numerov method
W = @(x) -x;

% control parameters
xmax = 5;
xmin = -15;

% Integrating from xmax to 0
N = 200;
h = (xmin-xmax)/N;
x = linspace(xmax,xmin,N+1);
y = zeros(1,N+1);
w = zeros(1,N+1);

% initial conditions
y(1) = 0;
w(1) = W(x(1));

% we guess next value
y(2) = y(1)+0.1;
w(2) = W(x(2));

% shoot left by the Numerov method
for n=2:N
  w(n+1) = W(x(n+1));
  y(n+1) = 2*(1-5*h^2*w(n)/12)*y(n) - (1+h^2*w(n-1)/12)*y(n-1); 
  y(n+1) = y(n+1)/(1+h^2*w(n+1)/12);
end

% normalization
N0 = int32(-xmax/h)+1; % find the location of x=0
y = y/y(N0);


p=plot(x,y,x,airy(x)/airy(0));
xlabel('x')
ylabel(texlabel('psi(x)'))
set(p(1),'linewidth',2,'color','red');
set(p(2),'linewidth',1,'color','blue');
legend('Numerov','MATLAB');
legend('location','southeast');
hold on
q=plot([-15, 5],[0,0],'--',[0,0],[-1.5,2],'--');
set(q,'color','black');
axis([-15 5 -1.5 2]);
hold off



