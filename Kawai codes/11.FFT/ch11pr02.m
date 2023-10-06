%**************************************************************************
%*     Example 11.2                                                       *
%*     filename: ch11pr02.m                                               *
%*     program listing number: 11.2                                       *
%*                                                                        *
%*     This program calculates the fourier transform of the second order  *
%*     derivative.                                                        *
%*                                                                        *
%*     Uses MATLAB function ifft() and fft()                              *
%*                                                                        *
%*     Programed by Ryoichi Kawai for Computational Physics Course.       *
%*     Last modification:  02/18/2017.                                    *
%**************************************************************************
clear all;

% Number of grid points
N=1024;

% Setting position domain
X=50.;   % size of position domain
dx=X/N;  % position interval
xmin=-X/2.;
xmax=xmin+dx*(N-1);
x0=linspace(xmin,xmax,N);    % position space in natural order
x1=linspace(0.0,dx*(N-1),N); % position space in fft order 

% setting momentum domain
K=2.0*pi*N/X;    % size of the momentum domain
dk=K/N;          % resolution in momentum space
kmin=-K/2.0;
kmax=kmin+dk*(N-1);
k0=linspace(kmin,kmax,N);     % momentum space in natural order
k1=fftshift(k0);              % momentum space in shifted order

f0=exp(-x0.^2)/sqrt(pi); % function in normal order
f1=ifftshift(f0);          % function in fft order

% FFT
F1=ifft(f1);        % prefactor not need this time
F1=-F1.*(k1.^2);    % because we perform both forward
g1=fft(F1);         % and inverse transformation.
g0=fftshift(g1) ;    % output in natural order

g2=exp(-x0.^2)/sqrt(pi).*(4*x0.^2-2); % analytic FT

p=plot(x0,real(g0),'-o',x0,g2);
set(p(1),'linewidth',2,'color','red')
set(p(2),'color','blue')
axis([-6 6 -1.4 0.8])
xlabel('x','fontsize',14)
ylabel('g(x)','fontsize',14)
legend('FFT','exact')
legend('location','southeast')