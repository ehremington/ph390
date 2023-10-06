%**************************************************************************
%*     Section 5.5.5                                                      *
%*     filename: ch05pr11.m                                               *
%*     program listing number: 5.11                                       *
%*                                                                        *
%*     This program finds the trajectory of a pendulum using              *
%*     Euler and Verlet methods.    Euler method shows its numerical      *
%*    instability and the trajectory diverges.                            *
%*                                                                        *
%*     Programed by Ryoichi Kawai for Computational Physics Course        *
%*     Revised on 01/07/2014.                                             *
%**************************************************************************


% system parameters
mass=1; L=1; g=9.8; Omega=sqrt(g/L); I=mass*L^2;

% initial conditions
theta1(1)=0.5; omega1(1)=0;t(1)=0;
theta2(1)=0.5; omega2(1)=0;
E1(1) = - mass*g*L*cos(theta1(1));
E2(1) = - mass*g*L*cos(theta2(1));

% control parameters
tmax=50; N=5000; h=tmax/N;

% Euler method
for i=1:N-1
    t(i+1)=t(1)+i*h;
    omega1(i+1) = omega1(i) - Omega^2 * sin(theta1(i)) * h;
    theta1(i+1) = theta1(i) + omega1(i)*h;
    E1(i+1) = I/2 * omega1(i+1)^2 - mass*g*L*cos(theta1(i+1));
end

% Verlet method
theta2(2) = theta2(1) + omega2(1)*h- Omega^2 * sin(theta2(1))*h^2/2;
for i=2:N
    theta2(i+1)=2*theta2(i)-theta2(i-1)-Omega^2 * sin(theta2(i))*h^2;
    omega2(i) = (theta2(i+1)-theta2(i-1))/(2*h);
    E2(i)=I/2 * omega2(i)^2 - mass*g*L*cos(theta2(i));
end

subplot(1,2,1);
q=plot(t(1:N),theta1(1:N),t(1:N),theta2(1:N));
xlabel('t');
ylabel(texlabel('theta'));
set(q(1),'Color','blue','Linewidth',2);
set(q(2),'Color','red','Linewidth',2);
legend(q,{'Euler','Verlet'});
legend(q,'Location','SouthWest');

subplot(1,2,2);
p=plot(t(1:N),E1(1:N),t(1:N),E2(1:N));
xlabel('t');
ylabel('Energy');
set(p(1),'Color','blue','Linewidth',2);
set(p(2),'Color','red','Linewidth',2);
legend(p,{'Euler','Verlet'});
legend(p,'Location','NorthWest');