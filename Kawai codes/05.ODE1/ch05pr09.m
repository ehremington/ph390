%**************************************************************************
%*     Section 5.5.3                                                      *
%*     filename: ch05pr09.m                                               *
%*     program listing number: 5.9                                        *
%*                                                                        *
%*     This program solves the synchronization of two phase oscillators   *
%*     using 2nd-order Runge-Kutta methods.                               *
%*                                                                        *
%*     Programed by Ryoichi Kawai for Computational Physics Course        *
%*     Revised on 01/07/2014.                                             *
%**************************************************************************
clear all;

omega1=1.0;  omega2=2.2;

% Control parameters
tmax=20; N=2000; h=tmax/N;

% initial conditions
theta1(1)=pi; theta2(1)=0; t(1)=0;

% 2nd-order Runge-Kutta method
for n=1:N-1
    t(n+1) = t(1)+n*h;
    k1 = omega1 + sin(theta2(n)-theta1(n));
    l1 = omega2 - sin(theta2(n)-theta1(n));
    mid1 = theta1(n)+k1*h/2;
    mid2 = theta2(n)+l1*h/2;
    k2 = omega1 + sin(mid2-mid1);
    l2 = omega2 - sin(mid2-mid1);
    theta1(n+1)=theta1(n)+k2*h;
    theta2(n+1)=theta2(n)+l2*h;
end

% plot the trajectories of oscillators
subplot(1,2,1); 
p=plot(t,sin(theta1),t,sin(theta2));
xlabel('t');
ylabel(texlabel('sin theta'));
set(p(1),'Color','blue','Linewidth',2);
set(p(2),'Color','red','Linewidth',2);
legend(p,{texlabel('theta_1'),texlabel('theta_2')});
legend(p,'Location','SouthEast');

% plot the phase difference
subplot(1,2,2);
q=plot(t,theta1-theta2);
xlabel('t');
ylabel(texlabel('theta_1-theta_2'));
set(q,'Linewidth',2);
