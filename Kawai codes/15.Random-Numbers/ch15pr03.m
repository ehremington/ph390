%**************************************************************************
%*     Example 15.3                                                       *
%*     filename: ch15pr03.m                                               *
%*     program listing number: 15.3                                       *
%*                                                                        *
%*     This program generates Gaussian distributed random numbers         *
%*     using the Box-Muller method.                                       *
%*                                                                        *
%*     Programed by Ryoichi Kawai for Computational Physics Course.       *
%*     Last modification:  02/25/2017.                                    *
%**************************************************************************
clear all

N=10000000;
x = rand(N,2);
y(:,1) = sqrt(-2*log(x(:,1))).*cos(2*pi*x(:,2));
y(:,2) = sqrt(-2*log(x(:,1))).*sin(2*pi*x(:,2));
ymin=min(y(:));
ymax=max(y(:));
fprintf('the largest deviation=%f, %f\n',ymin,ymax)

h=histogram(y,201,'Normalization','pdf');
hold on
% true normal distribution
K=201;
xmin=floor(ymin);
xmax=ceil(ymax);
X=linspace(xmin,xmax,K);
F=exp(-X.^2/2)/sqrt(2*pi);
p=plot(X,F);
set(p,'color','red','linewidth',2)
xlabel('x','fontsize',14)
ylabel('\rho(x)','fontsize',14)
legend('Box-Muller','Exact')
axis([xmin xmax 0.0 0.5])
hold off
