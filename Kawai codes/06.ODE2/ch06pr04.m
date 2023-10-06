%**************************************************************************
%*     Section 6.3.2                                                      *
%*     filename: ch06pr04.m                                               *
%*     program listing number: 6.4-1                                      *
%*                                                                        *
%*     This program solves one-dimensional heat equation and finds        *
%*     temperature profile and heat energy transaction using              *
%*     Numerov integration.                                               *
%*     Use function: numerov_heat(TL,delta,L)                             *
%*                                                                        *
%*     Programed by Ryoichi Kawai for Computational Physics Course.       *
%*     Last modification:  10/13/2013.                                    *
%**************************************************************************

clear all;

% system parameters
global mu
mu = -10;
TL=10;  TR=0; L=1;
% tolerance
tol=1e-9;
% control variable 
found = false;

n=1;
% first guess of delta
y1(n) = 1;
% get the u(x)
[x, y] = numerov_heat(TL,y1(n),L);
N = size(y,2);   % check how many grid points are used.
y2(n)=(y(N)-TR)^2;
if abs(y2(n)) < tol
    found = true;
end

if not(found)
    n=n+1;
    % second guess of delta
    y1(n) = y1(n-1)+0.01;
    % get u(x)
    [x, y] = numerov_heat(TL,y1(n),L);
    y2(n)=(y(N)-TR)^2;
    if abs(y2(n)) < tol
         found = true;
    end
end

% secant iteration
while not(found)
    % guess delta by secant method
    y1(n+1) = y1(n) - (y1(n)-y1(n-1))/(y2(n)-y2(n-1))*y2(n);
    % derivative of phi(x) at the end point.
    [x, y] = numerov_heat(TL,y1(n+1),L);
    y2(n+1)=(y(N)-TR)^2;
    if abs(y2(n+1)) < tol
        found = true;
    end
    n=n+1;
end

% Energy conservation
Q_in = -(y(2)-y(1))/(x(2)-x(1));
Q_out= +(y(n)-y(n-1))/(x(n)-x(n-1));
Q_diss = mu*sum(y(1:2:n-2)+4*y(2:2:n-1)+y(3:2:n))*(x(2)-x(1))/3;
fprintf('Q_in=%.6f,  Q_out=%.6f,  Q_diss=%.6f,  Q_net=%.6f\n',...
    Q_in, Q_out, Q_diss, Q_in+Q_out+Q_diss);

% plot heat source
subplot(1,2,1)
p=plot(x,y);
set(p,'Linewidth',2,'color','red')
xlabel('s','fontsize',14)
ylabel(texlabel('u(s)'),'fontsize',14)
hold on

p2=plot([0,L],[0,0],'--');
set(p2,'color','black')
hold off

subplot(1,2,2)
% plot the error after each iteration.
q=semilogy([1:n],abs(y2),'-o');
xlabel('iteration','fontsize',14)
ylabel(texlabel('T_R'),'fontsize',14)
set(q,'linewidth',2)