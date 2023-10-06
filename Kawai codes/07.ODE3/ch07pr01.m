%**************************************************************************
%*     Example  7.1                                                       *
%*     filename: ch07pr01.m                                               *
%*     program listing number: 7.1                                        *
%*                                                                        *
%*     This program finds the first two eigen modes of standing wave in a *
%*     string using the shooting method (Numerov and secant methods).     *
%*                                                                        *
%*     Programmed by Ryoichi Kawai for Computational Physics Course.      *
%*     Last modification:  02/04/2017.                                    *
%**************************************************************************
clear all;

% setting the grid
xmin=0.0;
xmax=1.0;
N=200;
h=(xmax-xmin)/N;
x = linspace(xmin,xmax,N+1);
v=zeros(1,N+1);

% control parameters
tol=1e-8;
delta=1;
found=false;
kmax=100;


lambda(1) = 1;
k=1;

while not(found) & k<kmax

    % integrate ODE by Numerov method
    w = -lambda(k);
    v(1)=0;
    v(2)=delta;
    for n=2:N-1
        v(n+1) = 2*(1-5*h^2*w/12)*v(n) - (1+h^2*w/12)*v(n-1);
        v(n+1) = v(n+1)/(1+h^2*w/12);
    end

    % error in the boundary condition
    err(k) = v(N);

    if abs(err(k)) < tol
        found = true;
    else
        % secant method to guess next lambda
        if k == 1
            lambda(k+1) = lambda(k)-0.1;
        else
            lambda(k+1) = lambda(k) ...
                -(lambda(k)-lambda(k-1))/(err(k)-err(k-1))*err(k);
        end
        k=k+1;
    end
end
s=max(abs(v));
v=v/s;
fprintf('lambda = %.7f, exact=%.7f\n',lambda(k),-pi^2);


subplot(1,2,1)
p=plot(x,v);
set(p,'linewidth',2)
xlabel('x','fontsize',14);
ylabel('v(x)','fontsize',14);
axis([ 0 1 -1.1 1.1]);
hold on

subplot(1,2,2)
r=plot([1:k],lambda(1:k),'-o');
set(r,'linewidth',2)
xlabel('iteration','fontsize',14)
ylabel(texlabel('lambda'),'fontsize',14)
axis([0 k -(2*pi)^2*1.1 10] )
hold on
r=plot([0,k],[-pi^2, -pi^2],'--');
set(r,'color','black')
hold on
r=plot([0,k],[-(2*pi)^2, -(2*pi)^2],'--');
set(r,'color','black')
hold on

found = false;
lambda(1) = -30;
k = 1;

while not(found)
    w = -lambda(k);
    v(1)=0;
    v(2)=delta;
    for n=2:N-1
        v(n+1) = 2*(1-5*h^2*w/12)*v(n) - (1+h^2*w/12)*v(n-1);
        v(n+1) = v(n+1)/(1+h^2*w/12);
    end
    err(k) = v(N);
    if abs(err(k)) < tol
        found = true;
    else
        if k == 1
            lambda(k+1) = lambda(k)-0.1;
        else
            lambda(k+1) = lambda(k) ...
                -(lambda(k)-lambda(k-1))/(err(k)-err(k-1))*err(k);
        end
        k=k+1;
    end
end
s = max(v);
v = v/s;
fprintf('lambda = %.7f, exact=%.7f\n',lambda(k),-(2*pi)^2);

subplot(1,2,1)
p=plot(x,v);
set(p,'linewidth',2,'color','red')
legend(texlabel('lambda_0=1'),texlabel('lambda_0=-30'));
legend('location','southwest');
r=plot([0,1],[0,0]);
set(r,'color','black')
hold off
subplot(1,2,2)
r=plot([1:k],lambda(1:k),'-o');
set(r,'linewidth',2,'color','red')
hold off
