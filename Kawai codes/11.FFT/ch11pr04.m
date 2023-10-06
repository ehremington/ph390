%**************************************************************************
%*     Section 11.4.4                                                     *
%*     filename: ch11pr04.m                                               *
%*     program listing number: 11.4                                       *
%*                                                                        *
%*     This program calculates the probability distribution of a harmonic *
%*     oscillator in the momentum spapce.
%*                                                                        *
%*     Uses MATLAB function ifft()                                        *
%*                                                                        *
%*     Programed by Ryoichi Kawai for Computational Physics Course.       *
%*     Last modification:  02/08/2015.                                    *
%**************************************************************************
clear all

% Number of grid points
N=1024;

% Setting position domain
X=100.;   % size of position domain
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
k1=linspace(0.0,dk*(N-1),N);  % momentum space in shifted order

% function in natural order
psix0=exp(-x0.^2/2)/sqrt(sqrt(pi));
psix1=2*x0.*exp(-x0.^2/2)/sqrt(2*sqrt(pi));
psix2=(2*x0.^2-1).*exp(-x0.^2/2)/sqrt(2*sqrt(pi));

% Ground state (n=0)
f=ifftshift(psix0);       % input in fft order
g=ifft(f)*X;              % do not forget to multiply X
psik0=ifftshift(g);
% 1st excited state (n=1)
f=ifftshift(psix1);       % input in fft order
g=ifft(f)*X;              % do not forget to multiply X
psik1=ifftshift(g);
%  2nd excited state (n=2)
f=ifftshift(psix2);       % input in fft order
g=ifft(f)*X;              % do not forget to multiply X
psik2=ifftshift(g);

subplot(3,1,1)
p1=plot(k0(N/4:3*N/4),abs(psik0(N/4:3*N/4)).^2);
ylabel('$|\psi_0(k)|^2$')
subplot(3,1,2)
p2=plot(k0(N/4:3*N/4),abs(psik1(N/4:3*N/4)).^2);
ylabel('$|\psi_1(k)|^2$')
subplot(3,1,3)
p3=plot(k0(N/4:3*N/4),abs(psik2(N/4:3*N/4)).^2);
ylabel('$|\psi_2(k)|^2$')
xlabel('k')
hold off


