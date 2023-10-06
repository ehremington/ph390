%**************************************************************************
%*     Section 12.3.1                                                     *
%*     filename: ch12pr07.m                                               *
%*     program listing number: 12.7                                       *
%*                                                                        *
%*     This program finds the activation energy of a reaction from a data *
%*     set using the linear regression.                                   *
%*                                                                        *
%*     Programed by Ryoichi Kawai for Computational Physics Course.       *
%*     Last modification:  02/25/2017.                                    *
%**************************************************************************
clear all;

k=8.617e-5; % Boltzmann constant [eV/K]

% Generate a new dataset
%E=0.025;
%A=2.0;
%for i=1:11
%    T(i) = 180+20*i;
%    f(i)=2.0*exp(-E/(k*T(i)))*(1+random('unif',-0.05,+0.05));
%end

T=[200, 220, 240, 260, 280, 300, 320, 340, 360, 380, 400];
f=[0.471,0.515,0.576,0.639,0.734,0.742,0.833,0.830,0.932,0.918,0.939];
N=size(f,2);

for i=1:N
    x(i)=1/(k*T(N+1-i));
    y(i)=log(f(N+1-i));
end

F=sum(y);
X=sum(x);
X2=sum(x.^2);
XF=sum(x.*y);
b=(F*X2-X*XF)/(N*X2-X^2);
a=(N*XF-X*F)/(N*X2-X^2);
g=a*x+b;
A=exp(b);
z = A*exp(a./(k*T));

fprintf('Activation Energy = %6.2d eV\n',abs(a))

subplot(1,2,1)
p=plot(x,g);
set(p,'linewidth',2,'color','red')
hold on
r=plot(x,y,'o');
set(r,'linewidth',2,'color','black')
xlabel(texlabel('beta'),'fontsize',14)
ylabel('log k','fontsize',14)
hold off
subplot(1,2,2)
p=plot(T,f,'o');
set(p,'linewidth',2,'color','black')
hold on
r=plot(T,z);
set(r,'linewidth',2,'color','red')
xlabel('T','fontsize',14)
ylabel('k','fontsize',14)
hold off
