%**************************************************************************
%*     Example  6.2                                                       *
%*     filename: numerov_poisson.m                                        *
%*     program listing number: 6.2-2                                      *
%*     Called by ch06pr02.m                                               *
%*                                                                        *
%*     This function integrates a one-dimensional Poisson equation for    *
%*     given initial velocity and the final time.                         *
%*        Input:  y1 = y(h)                                               *
%*                 L = boundary                                           *
%*        Output:  y = y(L)                                               *
%*                                                                        *
%*     Programed by Ryoichi Kawai for Computational Physics Course.       *
%*     Last modification:  10/13/2013.                                    *
%**************************************************************************
function [x, y]=numerov_poisson(y1,L,N)

% define S(x) in Numerov method
S=@(x) -2*x*exp(-x^2);

% control parameters
h=L/N;
x=linspace(0,L,N+1);
y=zeros(1,N+1);

% initial conditions
% due to symmetry phi(0)=0
y(1)=0;
s(1)=S(x(1));

% we guess phi(h)=phi_1
y(2)=y1;
s(2)=S(x(2));

% shoot out to x=L by the Numerov method
for n=2:N
  s(n+1)=S(x(n+1));
  y(n+1) = 2*y(n)-y(n-1)+(s(n+1)+10*s(n)+s(n-1))*h^2/12;
end

end


