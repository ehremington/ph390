%**************************************************************************
%*     Section 6.3.2                                                      *
%*     filename: numerov_heat.m                                           *
%*     program listing number: 6.3-2                                      *
%*     Called by ch06pr04.m                                               *
%*                                                                        *
%*     This function integrates a one-dimensional heat equation.          *
%*        Input:  TL = temperature at the left end of the rod             *
%*             delta = decrease of the temperature at next point          *
%*                 L = length of the rod                                  *
%*        Output:  y = position y(:,1) and temperature profile y(:,2)     *
%*                                                                        *
%*     Programed by Ryoichi Kawai for Computational Physics Course.       *
%*     Last modification:  10/13/2013.                                    *
%**************************************************************************
function [x, y]=numerov_heat(TL,delta,L)

global mu

% control parameters
N=10000; h=L/N;
x=linspace(0,L,N+1);
y=zeros(1,N+1);

% define w(x) in Numerov method
w=mu;

% initial conditions
y(1)=TL;

% we guess u(h)
y(2)=y(1)-delta;

% shoot out to x=L by the Numerov method
for n=2:N
  y(n+1) = 2*(1-5*h^2*w/12)*y(n) - (1+h^2*w/12)*y(n-1);
  y(n+1) = y(n+1)/(1+h^2*w/12);
end

end


