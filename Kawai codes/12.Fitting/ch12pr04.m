%**************************************************************************
%*     Example 12.4                                                       *
%*     filename: ch12pr04.m                                               *
%*     program listing number: 12.4                                       *
%*                                                                        *
%*     This program interpolates 11-point data with the Lagrange          *
%*     polynomial method.                                                 *
%*                                                                        *
%*     Programed by Ryoichi Kawai for Computational Physics Course.       *
%*     Last modification:  02/25/2017.                                    *
%**************************************************************************
clear all

% data to be fitted
F=[0.0000, 0.6889, 0.6095, 0.0774, -0.3401, -0.3528,...
   -0.0842, 0.1620, 0.1997, 0.0681, -0.0736];
X=[0.0,1.0,2.0,3.0,4.0,5.0,6.0,7.0,8.0,9.0,10.];
N=size(F,2);

M=101;
z=linspace(0,X(N),M);
for j=1:M
    L=1;
    y(j)=0;
    for n=1:N
        L=1; % Lagrange basis polynomial
        for m=1:N
            if n~=m
                L=L*(z(j)-X(m))/(X(n)-X(m));
            end
        end
        y(j)=y(j)+L*F(n);
    end
end

v=sin(z).*exp(-0.2*z);

p=plot(X,F,'o',z,y,z,v,'--');
set(p,'linewidth',2);
xlabel('x','fontsize',14);
ylabel('f(x)','fontsize',14);
legend('Data','Lagrange polynomial','Original');
legend('location','northeast');





