%**************************************************************************
%*     Section 5.5.2                                                      *
%*     filename: ch05pr08.m                                               *
%*     program listing number: 5.8                                        *
%*                                                                        *
%*     This program solves the Maxwell-Bloch model of laser dynamics      *
%*     using Runge-Kutta 4th order methods.                               *
%*                                                                        *
%*     Programed by Ryoichi Kawai for Computational Physics Course        *
%*     Revised on 01/07/2014.                                             *
%**************************************************************************
clear all;

% system parametes uncomment ONLY the desired set
%Type A
%gamma1=0.1; gamma2=2; gamma3=3; c1=0.25; c2=0.2; c3=1;
%Type B
%gamma1=0.1; gamma2=10; gamma3=0.25; c1=1; c2=0.5; c3=1;
%Type C 
gamma1=1; gamma2=0.1; gamma3=0.25; c1=1; c2=0.1; c3=1;

lambda=input('Enter a value for lambda = ');

% Control parameters
tmax=500; N=5000; h=tmax/N;

% initial conditions
E(1)=1.0; P(1)=1.0; D(1)=1.0; t(1)=0;

% 2nd-order Runge-Kutta method
for n=1:N-1
    t(n+1)=t(1)+n*h;
    FE_n=-gamma1*E(n)+c1*P(n);
    FP_n=-gamma2*P(n)+c2*E(n)*D(n);
    FD_n=-gamma3*(D(n)-lambda)-c3*E(n)*P(n);
    E_mid = E(n)+FE_n*h/2;
    P_mid = P(n)+FP_n*h/2;
    D_mid = D(n)+FD_n*h/2;
    FE_mid=-gamma1*E_mid+c1*P_mid;
    FP_mid=-gamma2*P_mid+c2*E_mid*D_mid;
    FD_mid=-gamma3*(D_mid-lambda)-c3*E_mid*P_mid;
    E(n+1)=E(n)+FE_mid*h;
    P(n+1)=P(n)+FP_mid*h;
    D(n+1)=D(n)+FD_mid*h;
end

% plot the dynamics of E
subplot(1,2,1);
plot(t,E);
xlabel('t');
ylabel('E(t)');

% plot 3D phase trajectory
subplot(1,2,2);
plot3(E,D,P);
xlabel('E');
ylabel('D');
zlabel('P');
grid on
