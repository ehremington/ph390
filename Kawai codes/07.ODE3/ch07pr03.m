%**************************************************************************
%*     Section 7.2.2                                                      *
%*     filename: ch07pr03.m                                               *
%*     program listing number: 7.3                                        *
%*                                                                        *
%*     This program finds an eigenvalue and eigenfunction of a quantum    *
%*     bouncing ball within a given bracket using the shooting            *
%*     method (Numerov and secant methods).                               *
%*                                                                        *
%*     Programmed by Ryoichi Kawai for Computational Physics Course.      *
%*     Last modification:  02/04/2017.                                    *
%**************************************************************************
clear all;

% setting the grid
L=10; N=500; h=L/N;
x=linspace(0.0,L,N+1);
psi=zeros(N+1);

% control parameter
tol=1e-6;
delta = 0.01;
found = false;
kmax=100;


eigval(1) = input('Initial Guess =');  % initial guess

% secant iteration
k=1;
while not(found)
    % Numerov integration
    w = eigval(k)-x;
    psi(N+1)=0;
    psi(N)=delta;

    for n=N:-1:2
        psi(n-1) = 2*(1-5*h^2*w(n)/12)*psi(n) - (1+h^2*w(n+1)/12)*psi(n+1);
        psi(n-1) = psi(n-1)/(1+h^2*w(n-1)/12);
    end

    err(k) = psi(1);  % last point should be zero
    if abs(err(k)) < tol
        found = true;
    else
        if k == 1
            eigval(k+1) = eigval(k)-0.1;  % second guess
        else
            eigval(k+1) = eigval(k) ... % suggestion by the secant method
                -(eigval(k)-eigval(k-1))/(err(k)-err(k-1))*err(k);
        end
        k=k+1;
    end
end

% normalize the solution at the maximum
psi=psi/max(psi);
fprintf('lambda = %.7f\n',eigval(k));

% plot eigenfunction
subplot(1,2,1)
p=plot(x,psi);
set(p,'linewidth',2)
hold on
q=plot([x(1),x(N)],[0,0]);
set(q,'color','black')
xlabel('x','fontsize',14);
ylabel('v(x)','fontsize',14);
axis([ 0 L -1.1 1.1]);

% plot eigenvalue
subplot(1,2,2)
r=plot([1:k],eigval(1:k),'-o');
set(r,'linewidth',2)
xlabel('iteration','fontsize',14)
ylabel('E','fontsize',14)
axis([0 k -pi^2*1.1 10] )
hold on
r=plot([0,k],[-pi^2, -pi^2],'--');
set(r,'color','black')
hold off
