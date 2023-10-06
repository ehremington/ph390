%**************************************************************************
%*     Section 5.5.4                                                      *
%*     filename: ch05pr10.m                                               *
%*     program listing number: 5.10                                       *
%*                                                                        *
%*     This program solves Newton equation for simple harmonic oscillator *
%*     using Runge-Kutta 4th order methods.  Then, it determines the      *
%*     period of the oscillation.                                         *
%*                                                                        *
%*     Programed by Ryoichi Kawai for Computational Physics Course        *
%*     Revised on 02/27/2017.                                             *
%**************************************************************************
clear all;

% system parameter
omega=1; mass=1; spring_k=mass*omega^2;E=1; 

% control parameters
h=0.01;
Nmax=10;

% initial conditions

x0=0.0;
N=0; 
t=0.0;

% the first step
x1=x0;
v1=sqrt(2*(E-spring_k*x1^2/2)/mass);
x2 = x1 + v1*h - omega^2*x1*h^2/2;

while N<Nmax

    % Verlet method
    x3=2*x2-x1 - omega^2 * x2 * h^2;
    v = (x3-x1)/(2*h);
    t=t+h;
    
    % Check if it returned to the tarting point
    if x3-x0>0 & x2-x0< 0 
        N=N+1;
    end
    
    x1=x2;
    x2=x3;
end

x=x1-x0;
% adjustment of the return time
Fn = -omega^2*x/mass;
if Fn>0 
    delta = -2*x/(v+sqrt( v^2-2*x*Fn));
else
    delta = (-v-sqrt( v^2-2*x*Fn))/Fn;
end

tau = t+delta;
period = tau/N;

fprintf('Period: Verlet = %7.6f,  Exact = %7.6f \n',period,2*pi);



