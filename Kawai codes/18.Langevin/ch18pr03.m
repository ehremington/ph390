%**************************************************************************
%*     Exercise 18.3                                                      *
%*     filename: ch18pr03.m                                               *
%*     program listing number: 18.3                                       *
%*                                                                        *
%*     This program calculates temporal autocorrelation of velocity of    *
%*     one-dimensional Brownian particles.                                *
%*                                                                        *
%*     Programed by Ryoichi Kawai for Computational Physics Course.       *
%*     Last modification:  03/26/2017.                                    *
%**************************************************************************
close all
clc

% system parameters
m=1;  % mass
gamma=0.1/m; % friction
T=2; % temperature
w=sqrt(2*gamma*T)/m;

% control parameters
N=2^20; % number of data points
dt=0.01;
ds=sqrt(dt)*w;

% initial conditions
x(1)=0; % initial position
v(1)=1; % initial velocity

% random number
g=randn(N,1);

% solve Langevin Eq.
for i=1:N-1
    % Euler step
    x(i+1)=x(i)+v(i)*dt;
    v(i+1)=v(i)-gamma*v(i)*dt+g(i)*ds;
    %Heun step
    x(i+1)=x(i)+(v(i)+v(i+1))*dt/2;
    v(i+1)=v(i)-gamma*(v(i)+v(i+1))*dt/2+g(i)*ds;
end

% velocity autocorrelation
z=fft(v);
y=abs(z).^2;
vc=ifft(y)/N;
vc=vc/vc(1);
t=[0:N-1]*dt;
vx=exp(-gamma*t);

subplot(1,2,1)
p=plot(t,v);
hold on
p=plot([1,t(N)],[0,0],'--');
set(p,'color','red')
xlabel('t','fontsize',14)
ylabel('v(t)','fontsize',14)
axis([t(1) t(N) -7 7])
hold off

subplot(1,2,2)
nc=floor(50/dt);
p=plot(t(1:nc),vx(1:nc));
set(p,'linewidth',2,'color','red')
hold on
p=plot(t(1:nc),vc(1:nc));
set(p,'linewidth',2)
legend('Exact','Simulation');
legend('location','northeast')
xlabel('\tau','fontsize',14)
ylabel('<v(\tau) v(0)>','fontsize',14)
axis([0 t(nc) -0.2 1.2])
p=plot([0,t(nc)],[0,0],'--');
set(p,'color','black')
hold off