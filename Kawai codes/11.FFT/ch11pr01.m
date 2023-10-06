%**************************************************************************
%*     Example 11.1                                                       *
%*     filename: ch11pr01.m                                               *
%*     program listing number: 11.1                                       *
%*                                                                        *
%*     This program calculates the fourier transform of Gaussian.         *
%*                                                                        *
%*     Uses MATLAB function ifft()                                        *
%*                                                                        *
%*     Programed by Ryoichi Kawai for Computational Physics Course.       *
%*     Last modification:  02/08/2015.                                    *
%**************************************************************************
clear all

% control parameters
N=1024;

% Setting time domain
T=50.;   % size of time domain
dt=T/N;  % time interval
tmin=-T/2.;
tmax=tmin+dt*(N-1);
t0=linspace(tmin,tmax,N);    % original time
t1=linspace(0.0,dt*(N-1),N); % fft ordered time 

% setting frequency domain
W=2.0*pi*N/T; % size of the frequency domain
dw=W/N;          % frequency resolution
wmin=-W/2.0;
wmax=wmin+dw*(N-1);
w0=linspace(wmin,wmax,N);     % frequency domain we want
w1=linspace(0.0,dw*(N-1),N);  % fft frequency domain

% function in normal order
f0=exp(-t0.^2)/sqrt(pi);

% function in swapped order
f1=ifftshift(f0);

% FFT
g1=ifft(f1)*T;  % do not forget to multiply T

% swap back to normal order
g0=ifftshift(g1);

% analytic FT
g2=exp((-w0.^2)/4);

subplot(2,2,1)
p=plot(t0,f0);
set(p,'linewidth',2);
hold on

axis([-T/2*1.05 T*1.05 0 1.0]);
legend('Original data')
legend('location','northwest')
xlabel('t','fontsize',14)
ylabel('f(t)','fontsize',14)
q=plot([-T/2,-T/2],[0,1],'--',[T/2,T/2],[0,1],'--');
set(q,'color','black')
hold off

subplot(2,2,3)
p=plot(t1,f1);
set(p,'linewidth',2);
axis([-T/2*1.05 T*1.05 0 1.0]);
legend('Swapped data')
legend('location','northwest')
xlabel('t','fontsize',14)
ylabel('f(t)','fontsize',14)
hold on
q=plot([0,0],[0,1],'--',[T,T],[0,1],'--');
set(q,'color','black')
hold off

subplot(2,2,2)
p=plot(w1,real(g1));
set(p,'linewidth',2');
axis([-W/2*1.05 W*1.05 0 1.3])
legend('FT Raw')
legend('location','northwest')
xlabel('$\omega$','fontsize',14)
ylabel('$\tilde{f}(\omega)$','fontsize',14)
hold on
q=plot([0,0],[0,1.3],'--',[W,W],[0,1.3],'--');
set(q,'color','black')
hold off

subplot(2,2,4)
p=plot(w0,g0);
set(p,'linewidth',2');
axis([-W/2*1.05 W*1.05 0 1.3])
legend('FT Swapped')
legend('location','northwest')
xlabel('$\omega$','fontsize',14)
ylabel('$\tilde{f}(\omega)$','fontsize',14)
hold on
q=plot([-W/2,-W/2],[0,1.3],'--',[W/2,W/2],[0,1.3],'--');
set(q,'color','black')
hold off

figure
p=plot(w0,g0,'-o',w0,g2);
axis([-5 5 0 1.3])
set(p(1),'linewidth',2')
set(p(2),'color','red')
legend('FFT','Exact')
xlabel('$\omega$','fontsize',16)
ylabel('$\tilde{f}(\omega)$','fontsize',16)