%**************************************************************************
%*     Section 18.2.1                                                     *
%*     filename: ch18pr05.m                                               *
%*     program listing number: 18.5                                       *
%*                                                                        *
%*     This program simulates a flashing ratchet.                         *
%*                                                                        *
%*     Programed by Ryoichi Kawai for Computational Physics Course.       *
%*     Last modification:  03/26/2017.                                    *
%**************************************************************************
close all
clear all

% system parameters
gamma=0.1; %friction
T=1; % temperature
L=1; % periodicity
k=2*pi/L; % wave number
U0=10; % potential strength
F0=U0*pi/L; % force strength
Fx=2;  %external load
h=10; % pontential on-off period
D=T/gamma; %noise strength

% control parameter
tau=100;  % total time
dt=0.005;  % time step
ds=sqrt(2.0*D*dt);  % Wiener time step
N=ceil(tau/dt);  % number of iterations
M=1000;

% prepare the arrays
sumx1=zeros(1,N);
sumx2=zeros(1,N);
t=linspace(0,dt*(N-1),N);
x=zeros(1,N);

figure(1)
for j=1:M
    u=0;
    potential=true;
    
    x(1)=0; % initial position
    g=randn(1,N);  % generate Gaussian white noise
    
    for i=1:N-1
        u=u+dt;
        if potential  % diffusion inside the potential
            f1 = -F0*(2.0*cos(k*x(i))-cos(2.0*k*x(i)))+Fx;
            x(i+1) = x(i)+g(i)*ds + f1*dt;
            f2 = -F0*(2.0*cos(k*x(i+1))-cos(2.0*k*x(i+1)))+Fx;     
            x(i+1) = x(i)+g(i)*ds + (f1+f2)*dt/2.0;
        else  % free diffusion
            x(i+1) = x(i)+g(i)*ds+Fx*dt;
        end    
        if u>h
            potential = not(potential);
            u=0;
        end
    end
    if j < 11  % show only 10 trajectories
        plot(t,x)
        drawnow
        hold on
    end
    sumx1 = sumx1+x;
    sumx2 = sumx2+x.^2;
    
end

% statstical analysis
xmean=sumx1./M;
sumx2=sumx2./M;
sigma2=sumx2-xmean.^2;

p=plot([0,t(N)],[0,0],'--');
set(p,'color','red')
p=plot(t,xmean);
set(p,'color','red','linewidth',2)
xlabel('t','fontsize',14)
ylabel('x(t)','fontsize',14)
hold off

figure(2);
p=plot(t,xmean,t,sigma2);
set(p,'linewidth',2)
xlabel('t','fontsize',14)
ylabel('<x>, \sigma^2','fontsize',14)
legend('<x>','\sigma^2')
legend('location','east')
hold on
p=plot([0,tau],[0,0],'--');
set(p,'color','black')
hold off
