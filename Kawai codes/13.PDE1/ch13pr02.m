%**************************************************************************
%*     Example 12.2                                                       *
%*     filename: ch12pr02.m                                               *
%*     program listing number: 12.2                                       *
%*                                                                        *
%*     This program calculates quantum tunneling by solving a Shrodinger  *
%*     equation.  The Crank-Nicolson method is used.                      *
%*                                                                        *
%*     Programed by Ryoichi Kawai for Computational Physics Course.       *
%*     Last modification:  02/15/2014.                                    *
%**************************************************************************

clear all
close all

% system parameters
E=0.9;
k=sqrt(2*E);
a=5.0;
dt=0.1;

% control parameters
L=100;
h=0.05*min(a,2*pi/k);
N=2*round(L/h)+1;

% initial condition
x0=-L+5*a;
for j=1:N;
    x(j)=(j-(N+1)/2)*h;
    psi(j,1)=exp(-(x(j)-x0)^2/(2*a^2))*exp(i*k*x(j));
end
c=sum(abs(psi).^2)*h;
psi=psi/sqrt(c);

% construct matrix
A=zeros(N,N);
A(1,1)=complex(1,dt/h^2/2)/2;
A(1,2)=-i*dt/h^2/8;
for n=2:N-1
    A(n,n)=A(1,1);
    A(n,n-1)=A(1,2);
    A(n,n+1)=A(1,2);
end
A(N,N)=A(1,1);
A(N,N-1)=A(1,2);

% add potential barrier
V=complex(0,dt/4);
for n=(N+1)/2:(N+1)/2+5/h;
    A(n,n)=A(n,n)+V;
end

% solve Schrodinger equation.
for I=1:1000    
    chi = A\psi;
    psi = chi - psi;
    rho = abs(psi).^2;
    p=plot(x,rho); 
    set(p,'linewidth',2);
    hold on
    r=plot([0,0],[0,0.12],[5,5],[0,0.12],[0,5],[0.12,0.12]);
    set(r,'color','black');
    axis([-L L 0 0.12]);
    xlabel('$x$','interpreter','latex','fontsize',16)
ylabel('$|\psi(x)|^2$','interpreter','latex','fontsize',16)
hold off; drawnow;
end
    
% check the normalization
c=sum(abs(psi).^2)*h;
fprintf('Final Norm=%.6f\n',c)

% compute transmission/reflection probability
T=sum(abs(psi((N+1)/2+int32(5/h):N-1).^2))*h;
R=sum(abs(psi(1:(N-1)/2).^2))*h;
fprintf('Transmission Probability=%.6f\n',T)
fprintf('Reflection   Probability=%.6f\n',R)

