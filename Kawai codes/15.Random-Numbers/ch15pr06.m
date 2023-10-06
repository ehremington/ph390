%**************************************************************************
%*     Section 15.5.2                                                     *
%*     filename: ch15pr06.m                                               *
%*     program listing number: 15.6                                       *
%*                                                                        *
%*     This program generates distribution of particles under gravity     *
%*     thermal diffusion.                                                 *
%*                                                                        *
%*     Programed by Ryoichi Kawai for Computational Physics Course.       *
%*     Last modification:  02/26/2017.                                    *
%**************************************************************************
clear all

N=1000; % number of particles

x=rand(N,1); %horizontal position = uniform
y=rand(N,1); % vertical position = exponential
z=-log(y);
zmax=max(z(:))+1;

subplot(1,2,1)
p=plot([-0.005 1.005],[-0.05,-0.05], [-0.005,-0.005],...
       [-0.05,zmax],[1.005,1.005],[-0.05,zmax]);
set(p,'color','black','linewidth',2)
hold on
p=plot(x,z,'.');
set(p,'linewidth',2)
axis([-0.1 +1.1 -0.5 zmax])
xlabel('x','fontsize',14)
ylabel('z','fontsize',14)
hold off

subplot(1,2,2)
h=histogram(z,int32(2*zmax),'Normalization','pdf');
hold on
Z=linspace(0.0,zmax,201);
P=exp(-Z);
q=plot(Z,P,'--');
set(q,'linewidth',2,'color','red')
xlabel('z','fontsize',14)
ylabel('probability density','fontsize',14)
hold off
