%**************************************************************************
%*     Example  6.2                                                       *
%*     filename: ch06pr02.m                                               *
%*     program listing number: 6.2-1                                      *
%*                                                                        *
%*     This program solves one-dimensional Poisson equation using         *
%*     Numerov integration and secant root finding methods.               *
%*     Use function: numerov_poisson(y,L)                                 *
%*                                                                        *
%*     Programed by Ryoichi Kawai for Computational Physics Course.       *
%*     Last modification:  10/13/2013.                                    *
%**************************************************************************
clear all;

% set the boundary conditions
L=10;
N=10000;

% tolerance
tol=1e-16';
% control variable 
found = false;

n=1;
% first guess of phi_1
y1(n) = 0.1;

% get the potential phi(x)
[x, y] = numerov_poisson(y1(n),L,N);

% derivative of phi(x) at the end point.
y2(n) = (y(N)-y(N-1))/(x(N)-x(N-1));
if abs(y2(n)) < tol
    found = true;
end

if not(found)
    n=n+1;
    % second guess of phi_1
    y1(n) = y1(n-1)+0.01;
    % get the potential phi(x)
    [x, y] = numerov_poisson(y1(n),L,N);
    % derivative of phi(x) at the end point.
    y2(n) = (y(N)-y(N-1))/(x(N)-x(N-1));
    if abs(y2(n)) < tol
         found = true;
    end
end

% secant iteration
while not(found)
    % guess phi_1 by secant method
    y1(n+1) = y1(n) - (y1(n)-y1(n-1))/(y2(n)-y2(n-1))*y2(n);
    % derivative of phi(x) at the end point.
    [x, y] = numerov_poisson(y1(n+1),L,N);
     % derivative of phi(x) at the end point.
    y2(n+1) = (y(N)-y(N-1))/(x(N)-x(N-1));
    if abs(y2(n+1)) < tol
         found = true;
    end
    n=n+1;
end

% construct the whole curve from x=-L to x=L.plt.
X(1:N) = -x(N:-1:1); X(N+1:2*N)=x(1:N);
Y(1:N) = -y(N:-1:1); Y(N+1:2*N)=y(1:N);

%plot charge density
subplot(1,2,1)
p1=plot(X,2.*X.*exp(-X.*X));
set(p1,'color','black')
hold on
% plot the numerical potential
p2=plot(X,Y);
set(p2,'Linewidth',2,'color','red')
xlabel('x','fontsize',14)
ylabel(texlabel('phi(x)'),'fontsize',14)
axis([-L L -1 1])
hold on
% plot the analytic potential
p3=plot(X,sqrt(pi)/2*erf(X));
set(p3,'color','blue')
legend('charge','numerical','exact')
legend('Location','southeast')
hold on
% plot the zero line
p4=plot([-L,L],[0,0],'--');
set(p4,'color','black')
hold off

subplot(1,2,2)
% plot the improvment of the first point.
q=semilogy([1:n],abs(y2),'-o');
xlabel('iteration')
ylabel(texlabel('phi''(L)'))
%axis([0 8 0.0001 1])
set(q,'linewidth',2)