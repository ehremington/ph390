
%**************************************************************************
%*     Example  4.1                                                       *
%*     filename: ch04pr01.m                                               *
%*     program listing number: 4.1                                        *
%*                                                                        *
%*     This program finds roots of a quadratic equation                   *
%*               a x^2 + x + 1/4 = 0                                      *
%*     for which the quadratic equation formula                           *
%*            x = (-1 + sqrt(1 - a))/(2a)                                 *
%*     fails for small value of a.                                        *
%*                                                                        *
%*     Programed by Ryoichi Kawai for Computational Physics Course.       *
%*     Last modification:  10/13/2013.                                    *
%**************************************************************************
clear all;

b=1;  % fixed parameters
c=0.25;

n=1;
a(1)=1;
while (a>eps())
  d = b^2-4*a(n)*c;
  x1(n) = (-b+sqrt(d))/(2*a(n));
  x2(n) = -2*c/(b+sqrt(d));
  fprintf('a= %22.16e, lhs=%22.16e, rhs=%22.16e \n',a(n),x1(n),x2(n));
  n=n+1;
  a(n)=a(n-1)/10;
end

% plot the results
p=loglog(a(1:n-1),abs(x2+0.25),'o',a(1:n-1),abs(x1+0.25),'s');
set(p(1),'Linewidth',2,'Color','red')
set(p(2),'Linewidth',2,'Color','blue')
xlabel('a','Fontsize',14)
ylabel(texlabel('|x - x_0|'),'Fontsize',14)  % x_0 = root at a=0
legend('lhs','rhs')
legend('Location','Southeast')

