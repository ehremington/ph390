%**************************************************************************
%*     Exercise 16.2                                                      *
%*     filename: ch16pr02.m                                               *
%*     program listing number: 16.2                                       *
%*                                                                        *
%*     This program simulates the one-dimensional persistent random walk. *
%*                                                                        *
%*     Programed by Ryoichi Kawai for Computational Physics Course.       *
%*     Last modification:  03/04/2017.                                    *
%**************************************************************************
clear all;

% parameters
p=0.75;  % persistent jump probability
M=100000; % number of particles
N=1000; % number of steps

% initial states
x=zeros(M,1);
mu=zeros(N,1);
sigma2=zeros(N,1);
mu(1)=0.;
sigma2(1)=0.;

d=randi(2,[M,1])*2-3; %unbiased initial jump
x=x+d;
mu(2)=sum(x)/M;  % mean
sigma2(2)=sum(x.^2)/M-mu(2)^2; %variance

for i=3:N
    r=rand(M,1);
    k=r>p;  % direction is reversed.
    d(k)=-d(k);
    x=x+sign(d);  % jump
    mu(i)=sum(x)/M;  % mean
    sigma2(i)=sum(x.^2)/M-mu(i)^2; %variance
end

subplot(1,2,1)
p=plot([1:N],mu,[1:N],sigma2);
set(p(1),'linewidth',2,'color','blue')
set(p(2),'linewidth',2,'color','red')
hold on
p=plot([1:N],[1:N],'--');
set(p,'color','black')
legend('\mu','\sigma^2','\sigma^2 (Normal RW)')
legend('location','northwest')
xlabel('steps')
ylabel('moments')
hold off

subplot(1,2,2)
h=histogram(x,51,'Normalization','pdf');
y=h.BinEdges;
g=1.0/sqrt(2*pi*N)*exp(-y.^2/(2.0*N));
hold on
r=plot(y,g);
set(r,'color','red')
xlabel('x')
ylabel('probability density')
legend('persistent RW','normal RW')
hold off
