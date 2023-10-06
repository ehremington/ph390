%**************************************************************************
%*     Section 18.2.2                                                     *
%*     filename: ch18pr06.m                                               *
%*     program listing number: 18.6                                       *
%*                                                                        *
%*     This program simulates a stochastic resonace.                      *
%*                                                                        *
%*     Programed by Ryoichi Kawai for Computational Physics Course.       *
%*     Last modification:  03/26/2017.                                    *
%**************************************************************************
close all
clear all

% system parameters
% To see stochastic resonance, use the parameter valesu
% A=1.0; D=1.2; U0=10; w=0.5;
%
A=1.0;
D=1.2;
U0=10;
w=0.5;

% constrol parameters
movie=false;
tau=100;
N=2^14;
dt=tau/N;
ds=sqrt(2.*D*dt);

x=zeros(1,N);
t=linspace(0.0,dt*(N-1),N);

L=101;
q=linspace(-2.0,2.0,L);
U1=U0*(q.^4/4-q.^2/2);

if movie
    figure(1)
    plot(q,U1-A*cos(w*t(1))*q);
    hold on
    y=U0*(x(1)^4/4-x(1)^2/2)+0.05-A*cos(w*t(1))*x(1);
    rectangle('Position',[x(1)-0.05,y-0.05,0.1,0.1],'Curvature',[1,1],'FaceColor',[0,0.75,0.75])
    axis equal
    axis([-2 2 -U0/2 2])
    drawnow;
    hold off
end

g=randn(1,N);

for i=1:N-1
    t(i+1)=i*dt;
    f1=U0*(-x(i)^3+x(i))+A*cos(w*t(i));
    x(i+1)=x(i) + f1*dt + g(i)*ds;
    f2=U0*(-x(i+1)^3+x(i+1))+A*cos(w*t(i+1));
    x(i+1)=x(i) + (f1+f2)*dt/2 + g(i)*ds;
    if movie 
        plot(q,U1-A*cos(w*t(i+1))*q);
        hold on
        y=U0*(x(i+1)^4/4-x(i+1)^2/2)+0.05-A*cos(w*t(i+1))*x(i+1);
        rectangle('Position',[x(i+1)-0.05,y-0.05,0.1,0.1],'Curvature',[1,1],'FaceColor',[0,0.75,0.75])
        axis equal
        axis([-2 2 -U0/2 2])
        drawnow;
        hold off
    end
end

figure(2)
plot(t,x)
hold on
p=plot(t,cos(w*t),'--');
set(p,'color','red')
axis([t(1),t(N),-2 2])
xlabel('t','fontsize',14')
ylabel('x(t)','fontsize',14)
hold off

figure(3)
z=fft(x)*tau/N;
q=linspace(0,N-1,N)*2*pi/tau;
p=plot(q,abs(z).^2);
set(p,'linewidth',2)
xlim([0 3])
ylabel('Power spectrum','fontsize',14)
xlabel('\omega','fontsize',14)
hold on
p=plot([0.5,0.5],[0,900],'--');
set(p,'color','black')
hold off
