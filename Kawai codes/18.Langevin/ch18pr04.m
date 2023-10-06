%**************************************************************************
%*     Example 18.4                                                       *
%*     filename: ch18pr04.m                                               *
%*     program listing number: 18.4                                       *
%*                                                                        *
%*     This program simulates the Orstein-Uhlenbeck process.              *
%*                                                                        *
%*     Programed by Ryoichi Kawai for Computational Physics Course.       *
%*     Last modification:  03/26/2017.                                    *
%**************************************************************************
clear all
clc

% system parameters
gamma=0.1; %friction
T=1.0; %temperature
k=1.0; %spring constant

% control parameter
tau=1000;
dt=0.01;
ds=sqrt(2*T/gamma*dt);
dk=k*dt;
N=ceil(tau/dt); % numner of time steps
M=1000; %number of realization

t=(0:N-1)*dt;

sumx1=zeros(N,1);
sumx2=zeros(N,1);

for j=1:M
    
    x(1)=0;  % initial position
    r=rand(N,2);                                 % Gaussian random numner
    g = sqrt(-2*log(r(:,1))).*cos(2*pi*r(:,2));  % by the Box-Muller

    for i=2:N
        t(i)=(i-1)*dt;
        x(i) = x(i-1)+ds*g(i) - dk*x(i-1);  % Eulrer step
        x(i) = x(i-1)+ds*g(i) - dk*(x(i)+x(i-1))/2;  % Huen step
    end
    sumx1(:) = sumx1(:)+x(:);
    sumx2(:) = sumx2(:)+x(:).^2;
    
end

figure(1);
plot(t,x)
hold on
p=plot([0,t(N)],[0,0],'--');
set(p,'color','red')
xlabel('t','fontsize',14)
ylabel('x(t)','fontsize',14)
hold off

figure(2);
xmean=sumx1./M;
sumx2=sumx2./M;
sigma2=sumx2-xmean.^2;
theory = T/gamma;
p=plot(t,xmean,t,sigma2,[t(1),t(N)],[theory,theory]);
set(p,'linewidth',2)
xlabel('t','fontsize',14)
ylabel('<x>, \sigma^2','fontsize',14)
legend('<x>','\sigma^2','theory')
legend('location','east')

figure(3)
[f,z]=hist(x,21);
dz=z(2)-z(1);
s=sum(f)*dz;
f=f/s;
y0=z(1);
y1=z(end);
y=[1:100]*(y1-y0)/100;
y=y+y0;
h=sqrt(gamma/(2*pi*T))*exp(-k*y(:).^2*gamma/(2*T));
p=plot(z,f,'s',y,h);
xm=max(abs(z(:)));
set(p(1),'MarkerFaceColor','black')
set(p(2),'color','red','linewidth',2)
legend('simulation','theory')
xlabel('x','fontsize',14)
ylabel('\rho(x)','fontsize',14)
axis([-xm xm 0 0.2])

