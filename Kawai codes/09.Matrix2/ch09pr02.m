%**************************************************************************
%*     Section  9.3.1                                                     *
%*     filename: ch09pr02.m                                               *
%*     program listing number: 9.2-1                                      *
%*                                                                        *
%*     This program finds a fixed point of Maxwell-Bloch equation         *
%*     by newton-raphson method and Gaussian elimination.                 *
%*                                                                        *
%*     Programed by Ryoichi Kawai for Computational Physics Course.       *
%*     Last modification:  02/10/2017.                                    *
%**************************************************************************
clear all;

% system parameters
g1=0.1; g2=2; g3=3;
k1=0.25; k2=0.2; k3=1;
lambda=5;

% control parameters
alpha = 0.5;
tol = 1e-4;
N=3;

%iniital guess
y(1,1:3)=[2;2;5];
i=1;
err = 1;

x(1:3)=y(i,:);
J=[[-g1,k1,0];[k2*x(3),-g2,k2*x(2)];[-k3*x(2),-k3*x(1),-g3]];
f=[-g1*x(1)+k1*x(2);-g2*x(2)+k2*x(1)*x(3);-g3*(x(3)-lambda)-k3*x(1)*x(2)];
while err > tol
    i=i+1;
    B=J\f;
    y(i,:)=y(i-1,:) - alpha*B';
    x=y(i,:);
    J=[[-g1,k1,0];[k2*x(3),-g2,k2*x(2)];[-k3*x(2),-k3*x(1),-g3]];
    f=[-g1*x(1)+k1*x(2);-g2*x(2)+k2*x(1)*x(3);-g3*(x(3)-lambda)-k3*x(1)*x(2)];
    err = f'*f;
end

fprintf('E=%.6f,  P=%.6f,  D=%.6f\n',x)

p=plot([1:i],y(:,1),'o-',[1:i],y(:,2),'o-',[1:i],y(:,3),'o-');
set(p,'linewidth',2)
xlabel('iteration','fontsize',14)
ylabel('lasing state','fontsize',14)
legend('E','P','D')
legend('location','northeast')