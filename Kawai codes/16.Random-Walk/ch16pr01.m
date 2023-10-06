%**************************************************************************
%*     Exercise 16.1                                                      *
%*     filename: ch16pr01.m                                               *
%*     program listing number: 16.1                                       *
%*                                                                        *
%*     This program simulates the one-dimensional descrete random walk.   *
%*                                                                        *
%*     Programed by Ryoichi Kawai for Computational Physics Course.       *
%*     Last modification:  03/03/2017.                                    *
%**************************************************************************
clear all

N=1000; % max time (number of steps)
M=100000; % number of particles

x=zeros(M,N); % reset the trajectories

% trajectory calculation
for i=1:N-1
    x(:,i+1)=x(:,i)+randi(2,[M,1])*2-3;
end

% stattistics
mu=sum(x,1)/M;  % mean
sigma=sqrt(sum(x.^2,1)/M);  %variance
t=[1:N];

subplot(1,2,1)
p=plot(t,mu,t,sigma,'--',t,-sigma,'--');
set(p(1),'color','red','linewidth',2)
set(p(2:3),'color','black','linewidth',2)
legend('\mu','+\sigma','-\sigma')
legend('location','northwest')
hold on
p=plot(t,x(1,:),t,x(2,:),t,x(3,:));
set(p(1),'color','blue');
set(p(2),'color','green');
set(p(3),'color',[0,0.75,0.75]);
axis([0 N -40 40])
xlabel('steps','fontsize',14)
ylabel('x','fontsize',14)
hold off

subplot(1,2,2)
rho=x(:,N);
h=histogram(rho,51,'Normalization','pdf');
hold on
X=h.BinEdges;
Y=1./sqrt(2*pi*N)*exp(-X.^2/(2.*N));
p=plot(X,Y);
set(p,'color','red','linewidth',2);
xlabel('x')
ylabel(texlabel('rho(x)'))
hold off