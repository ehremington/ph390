%**************************************************************************
%*     Exercise 18.2                                                      *
%*     filename: ch18pr02.m                                               *
%*     program listing number: 18.2                                       *
%*                                                                        *
%*     This program simulates the two-dimensional free continuous-time    *
%*     random walk (Wiener process) using the Heun method.                *
%*                                                                        *
%*     Programed by Ryoichi Kawai for Computational Physics Course.       *
%*     Last modification:  03/25/2017.                                    *
%**************************************************************************
close all
clc

% system parameters
gamma=0.1;  % frictional constant
T=1.0;      % temperature  
D=T/gamma;  % strength of noise

% control parameter
tau=100; % total time
dt=0.005; % time step
ds=sqrt(2.*D*dt); % Wiener step size
N=ceil(tau/dt); % number of steps
t=linspace(0.0,double(N*dt),N+1);
M=100; % number of samples

% reset the counters
sumx1=zeros(N+1,1);
sumx2=zeros(N+1,1);
sumy1=zeros(N+1,1);
sumy2=zeros(N+1,1);
x=zeros(N+1,1);
y=zeros(N+1,1);

% sampling loop begins
for j=1:M
    
    x(1)=0; % initial position
    y(1)=0;
    
    % Normally ditributed random numbers by the Box-Muller method
    g=randn(N,2);

    % integration of the Langevin equation
    for i=1:N
        x(i+1) = x(i)+g(i,1)*ds;
        y(i+1) = y(i)+g(i,2)*ds;
    end
    
    % record the data for statistical analysis
    sumx1 = sumx1+x;
    sumy1 = sumy1+y;
    sumx2 = sumx2+x.^2;
    sumy2 = sumy2+y.^2;
end

% plot a sample trajectory
figure(1)
plot(x,y)
hold on
plot(x(1),y(1),'.','color','Red')
axis equal;
xlabel('x','fontsize',14)
ylabel('y','fontsize',14)
hold off

% statistical analysis
figure(2);
xmean=sumx1/M;   % <x>
ymean=sumy1/M;   % <y>
sumx2=sumx2/M;   % <x^2>
sumy2=sumy2/M;   % <y^2>
xvar=sumx2-xmean.^2; % variance_x
yvar=sumy2-ymean.^2; % variance_y
theory=2*T/gamma*t;  % theoretical variance
p=plot(t,xmean,t,ymean,t,xvar,t,yvar,t,theory);
set(p,'linewidth',2)
xlabel('t','fontsize',14)
ylabel('<x>, \sigma^2','fontsize',14)
legend('<x>','<y>','\sigma_x^2','\sigma_y^2','theory')
legend('location','east')

