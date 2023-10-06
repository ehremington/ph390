%**************************************************************************
%*     Example 4.5                                                        *
%*     filename: ch04pr05.m                                               *
%*     program listing number: 4.5                                        *
%*                                                                        *
%*     This program calculates magnetization as a function of             *
%*     temperature using the mean field Ising model:                      *
%*         x = x tan(x/S)                                                 *
%*     where the variables are normalized as                              *
%*         x = m/m_0 and S=k*T/(C*m_0) .                                  *
%*                                                                        *
%*     Programed by Ryoichi Kawai for Computational Physics Course        *
%*     Revised on 01/07/2014.                                             *
%**************************************************************************
clear all;
N=101; % number of data points
tol=10^(-6);  % tolerance

dS = 2/(N-1); % step size in temperature

for i=1:N
    S=(i-1)*dS;  % temperature
    a=1/S;
    f=@(x) x - tanh(a*x);  % function
    df=@(x) 1 - a*sech(a*x)^2;  % derivative
    % Newton-Raphson method
    x=2;  % initial guess
    fx=f(x);
    while fx>tol
      x=x-fx/df(x);
      fx=f(x);
    end
    % store the magnetization and temperature  
    m(i)=x;
    t(i)=S;
end

p=plot(t,m,t,-m,[0,1],[0,0]);
xlabel('T','fontsize',14);
ylabel('m','fontsize',14);
set(p(1),'linewidth',2,'color','blue');
set(p(2),'linewidth',2,'color','blue');
set(p(3),'linewidth',2,'color','red');
axis([0 2 -1.1 1.1]);

