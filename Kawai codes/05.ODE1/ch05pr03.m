%**************************************************************************
%*     Example 5.3                                                        *
%*     filename: ch05pr03.m                                               *
%*     program listing number: 5.3                                        *
%*                                                                        *
%*     This program solves Newton equation for a falling object           *
%*     using the Runge-Kutta-Fehlberg method.                             *
%*        m = mass of the object                                          *
%*        g = acceleration due to gravity                                 *
%*        gamma = frictional coefficient                                  *
%*                                                                        *
%*     Programed by Ryoichi Kawai for Computational Physics Course        *
%*     Revised on 01/22/2017.                                             *
%**************************************************************************
clear all;

% system parameters
gamma=1.0;
g=9.8;
m=1.0;

tmax=10; % maximum time

% define the right hand side
f = @(t,v) -gamma*v - m*g;

% use RK45 method
[t,v]=ode45(f,[0,tmax],0.0);

% exact solution
v_ex = m*g/gamma * (exp(-gamma*t)-1);

subplot(1,2,1);
q=plot(t,v,'o',t,v_ex,'-');
xlabel('t','fontsize',14);
ylabel('v(t)','fontsize',14);
set(q(1),'Color','red','Linewidth',2);
set(q(2),'Color','black','Linewidth',2);
legend(q,'RK45','Exact')
legend(q,'Location','NorthEast');
hold on
plot(t,0,'ok')
hold off

subplot(1,2,2);
% Plot the absolute errors
p=semilogy(t,abs(v-v_ex));

% Plotting options
xlabel('t','fontsize',14);
ylabel('absolute error','fontsize',14);
set(p(1),'Color','blue','Linewidth',2);

