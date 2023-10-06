%**************************************************************************
%*     Example 12.1                                                       *
%*     filename: ch12pr01.m                                               *
%*     program listing number: 12.1                                       *
%*                                                                        *
%*     This program solves a diffusion equation using the forward time    *
%*     centered space method.                                             *
%*                                                                        *
%*     Programed by Ryoichi Kawai for Computational Physics Course.       *
%*     Last modification:  02/15/2014.                                    *
%**************************************************************************
close all;
clear all;
% parameters
D=0.1; % diffusion constant
N=201; % number of grids
N1=(N-1)/2;
dx=0.1; % spacial step
R = 0.1; % D*dt/dx^2
dt=dx^2/D *R;  % time step

x=(-N1:N1)*dx; % spatial coordinates

M=10; % number of sample point.
tmax=10000; % total time steps
MS=tmax/M;

rho0=zeros(N,1);  % initial density profile
rho0(N1+1)=1.0/dx;

% allocate arrays
rho1=zeros(N,1);
rho=zeros(N,M);
t=zeros(M,1);

k=0;
for j=1:tmax
    rho1(1)=(1-R)*rho0(1)+R*rho0(2);  % left boundary
    rho1(N)=(1-R)*rho0(N)+R*rho0(N-1);  % right boundary
    for i=2:N-1
        rho1(i)=rho0(i)*(1-2*R)+R*(rho0(i+1)+rho0(i-1));
    end
    rho0=rho1;
    if mod(j,MS)==0  % record the results
        k=k+1;
        t(k)=j*dt;
        rho(:,k)=rho0(:);
    end
end

figure(1)
surf(t,(-N1:N1),rho)
    
figure(2)
for i=1:M
    plot(x,rho(:,i))
    hold on
end
hold off

figure(3)
f=1/sqrt(2*pi)*1/sqrt(2*D*t(2))*exp(-x.^2/(4*D*t(2))); % exact
plot(x,rho(:,2),'o',x,f)
legend('FTCS method','Exact')