%**************************************************************************
%*     Section 16.4.3                                                     *
%*     filename: ch16pr06.m                                               *
%*     program listing number: 16.6                                       *
%*                                                                        *
%*     This program simulates the Parrondo game.                          *
%*                                                                        *
%*     Programed by Ryoichi Kawai for Computational Physics Course.       *
%*     Last modification:  03/04/2017.                                    *
%**************************************************************************
clear all
close all

% control parameters
e=0.005;
pA=1/2-e;
qA=1-pA;
pB=3/4-e;
qB=1-pB;
pC=1/10-e;
qC=1-pC;
N=50000;
M=100;

% Game A
x=zeros(N,1);
y=zeros(M+1,1);
for i=1:M
    r=rand(N,1);
    x(r<pA)=x(r<pA)+1;
    x(r>=pA)=x(r>=pA)-1;    
    y(i+1)=sum(x)/N;
end

p=plot([0:M],y);
set(p,'linewidth',2)
hold on

% Game B
x=zeros(N,1);
p=zeros(N,1);
for i=1:M
    r=rand(N,1);
    p(:)=pB;
    k=find(mod(x,3)==0);
    p(k)=pC;
    x(r<p)=x(r<p)+1;
    x(r>=p)=x(r>=p)-1;  
    y(i+1)=sum(x,1)/N;
end
p=plot([0:M],y);
set(p,'linewidth',2,'color',[0, 0.75,0.75])
hold on

% Game A and B
x=zeros(N,1);
p=zeros(N,1);
for i=1:M
    r=rand(N,1);
    if mod(i,4)<2
        p(:)=pA;
    else
        p(:)=pB;
        k=find(mod(x,3)==0);
        p(k)=pC;
    end
    x(r<p)=x(r<p)+1;
    x(r>=p)=x(r>=p)-1;  
    y(i+1)=sum(x,1)/N;
end
p=plot([0:M],y);
set(p,'linewidth',2,'color','red')
legend('Game A alone','Game B alone', 'Games A & B')
legend('location','northwest')
xlabel('# of plays','fontsize',14)
ylabel('Gain','fontsize',14)
hold off



    