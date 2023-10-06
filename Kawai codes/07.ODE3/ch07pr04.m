%**************************************************************************
%*     Section 7.2.3                                                      *
%*     filename: ch07pr04.m                                               *
%*     program listing number: 7.4                                        *
%*                                                                        *
%*     This program finds an eigenvalue and eigenfunction of a quantum    *
%*     particle in the Morse potential using the shooting.                *
%*     method (Numerov and secant methods).                               *
%*                                                                        *
%*     Programmed by Ryoichi Kawai for Computational Physics Course.      *
%*     Last modification:  10/13/2013.                                    *
%**************************************************************************
clear all
clc

% Use the atomic unit
mp=1836.4;
e=27.2114;
a0=0.529177;
% parameters for H2
D=4.75/e; R0=0.742/a0; a=1.44/R0; name='H_2'; m=mp; mu=m/2;

% parameters for I2
%D=1.56/e; R0=2.66/a0; a=4.94/R0; name='I_2'; m=127*mp; mu=m/2;

% zero point energy
omega=sqrt(2*D*a^2/mu);
DE=omega/5;
d=1/sqrt(mu*omega)*12;


% define the discrete coordinate
N=1000;
h=d/N;
x=linspace(-d/2,d/2,N+1);
psi=zeros(1,N+1);
phi=zeros(3,N+1);
% evaluate the potential
U=D*(exp(-2*a*x)-2*exp(-a*x));

E_M=-D;

for L=1:3
% Initial Bracketting
E_L = E_M+DE;
w = 2*mu*(E_L-U);
psi(1)=0;
psi(2)=0.1;
% shoot out to x=L by the Numerov method
for n=2:N
   psi(n+1) = 2*(1-5*h^2*w(n)/12)*psi(n) - (1+h^2*w(n-1)/12)*psi(n-1);
   psi(n+1) = psi(n+1)/(1+h^2*w(n+1)/12);
end
ERR_L=psi(N+1);
found=false;

while not(found)
    E_U = E_L+DE;
    if E_U>0
        error
    end
    w = 2*mu*(E_U-U);
    psi(1)=0;
    psi(2)=0.1;
    for n=2:N
        psi(n+1) = 2*(1-5*h^2*w(n)/12)*psi(n) - (1+h^2*w(n-1)/12)*psi(n-1);
        psi(n+1) = psi(n+1)/(1+h^2*w(n+1)/12);
    end
    ERR_U=psi(N+1);
    if ERR_L*ERR_U<0
        found=true;
    else
        E_L=E_U;
        ERR_L=ERR_U;
    end
end

found=false;
tol1=1e-8;
tol2=1e-15;

while not(found)
    E_M=(E_L+E_U)/2;
    w = 2*mu*(E_M-U);
    psi(1)=0;
    psi(2)=0.1;

    for n=2:N
        psi(n+1) = 2*(1-5*h^2*w(n)/12)*psi(n) - (1+h^2*w(n-1)/12)*psi(n-1);
        psi(n+1) = psi(n+1)/(1+h^2*w(n+1)/12);
    end
    ERR_M=psi(N+1);
    if abs(ERR_M) < tol1 || abs(E_L-E_M)<tol2
        found = true;
    end

    if ERR_M*ERR_L > 0
        E_L=E_M;
    else
        E_U=E_M;
    end

end
  eigval(L) = E_M;
  % normalize the wave function
  z=sum(psi.^2)*h;
  phi(L,:)=psi/sqrt(z);
  exact=-D*(1-a/sqrt(2*mu*D)*(L-1/2))^2;
  fprintf('n=%d, E=%.10f, Exact=%.10f\n',L-1,E_M,exact)
end


% Ploting potential and eigenvalues
subplot(1,2,1)
N0=ceil(N/5);
p=plot([x(N0),x(N)],[eigval(1),eigval(1)],...
       [x(N0),x(N)],[eigval(2),eigval(2)],...
       [x(N0),x(N)],[eigval(3),eigval(3)]);
legend(p,'n=0','n=1','n=2')
hold on
plot(x(N0:N),U(N0:N),'color','black')
hold off
xlabel(texlabel('xi'),'fontsize',14)
ylabel(texlabel('U(xi)'),'fontsize',14)

% Ploting wave functions
subplot(1,2,2)
plot(x,phi(1,:),x,phi(2,:),x,phi(3,:))
ymax=max(phi(:));
ymin=min(phi(:));
axis([-d/2 d/2 ymin*1.1 ymax*1.1])
xlabel(texlabel('xi'),'fontsize',14)
ylabel(texlabel('psi(xi)'),'fontsize',14)

