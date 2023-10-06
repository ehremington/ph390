%**************************************************************************
%*     Exercise 18.1                                                      *
%*     filename: ch18pr01.m                                               *
%*     program listing number: 18.1                                       *
%*                                                                        *
%*     This program simulates the one-dimensional free continuous-time    *
%*     random walk (Wiener process) using the Heun method.                *
%*                                                                        *
%*     Programed by Ryoichi Kawai for Computational Physics Course.       *
%*     Last modification:  03/25/2017.                                    *
%**************************************************************************
close all
clc

% system parameters
gamma=0.1; % friction coefficient
T=1.0; % temperature

% control parameter
tau=100;  % total duration
dt=0.005; % time step
D=T/gamma;  % diffusion constant
ds=sqrt(2.0*D*dt); % time step for Wiener process
N=int32(tau/dt);  % number of steps
M=1000; %number of samplings
t=linspace(0.0,double(N*dt),N+1);

figure(1)
theory=2.0*T/gamma*t; % theoretical variance
p=plot([t(1),t(N+1)],[0,0],'--');
set(p,'color','red')
p=plot(t,sqrt(theory),t,-sqrt(theory));
set(p,'color','red','linewidth',2)
xlabel('t','fontsize',14)
ylabel('x(t)','fontsize',14)
hold on

% reset to zero
sumx1=zeros(N+1,1);
sumx2=zeros(N+1,1);
x=zeros(N+1,1);

% loop over M samplings
for j=1:M
    % initial condition
    x(1)=0;
    % random force with Gaussian stribution
    g=randn(N,1);  

    % solving Langevin equation
    for i=1:N
        x(i+1) = x(i)+g(i)*ds;
    end
    if j<11
        plot(t,x)
        drawnow
    end
    
    % save the data for statistical analysis
    sumx1(:) = sumx1(:)+x(:);
    sumx2(:) = sumx2(:)+x(:).^2;
end

hold off

% statistical analysis
xmean=sumx1./M;   % <x>
sumx2=sumx2./M;   % <x^2>
xvar=sumx2-xmean.^2;  % variance

figure(2)
p=plot(t,xmean,t,xvar,t,theory)
set(p,'linewidth',2)
xlabel('t');
ylabel('<x>, \sigma^2','fontsize',14)
legend('\langle x \rangle','\sigma^2','exact')
legend('location','northwest')
hold off
