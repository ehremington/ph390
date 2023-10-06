%**************************************************************************
%*     Section 12.3.2                                                     *
%*     filename: ch12pr08.m                                               *
%*     program listing number: 12.8                                       *
%*                                                                        *
%*     This program finds the peak position and life-time broadening      *
%*     of atomic emmision spectrum.                                       *
%*                                                                        *
%*     Programed by Ryoichi Kawai for Computational Physics Course.       *
%*     Last modification:  02/25/2017.                                    *
%**************************************************************************
clear all;

% Generate a sample noisy data

x=[-2.01,-1.47,-0.97,-0.52,-0.04,0.52,0.99,1.53,2.03,2.51,2.96,3.47,4.02];
y=[0.28,0.57,0.62,0.68,1.26,1.29,1.57,1.11,0.91,0.94,0.65,0.80,0.31];
s=[0.10,0.11,0.17,0.06,0.15,0.11,0.15,0.10,0.11,0.14,0.16,0.18,0.15];
N=size(y,2);

%control parametrs
alpha=1e-2;
found=false;

%initial guess
M=3;
lambda(:,1)=[1,0,1];

% Gauss-Newton iteration
n=1;
% evaluate initial chi sqaure
for i=1:N
   F=lambda(1,n)/((x(i)-lambda(2,n))^2+lambda(3,n));
   b(i)=y(i)-F;
end
b=b./s;
chi2(n)=b*b';

while not(found)
    % construct Jacobian and vectors.
    for i=1:N
        F=lambda(1,n)/((x(i)-lambda(2,n))^2+lambda(3,n));
        J(i,1)=F/lambda(1,n);
        J(i,2)=2*lambda(1,n)*(x(i)-lambda(2,n))...
              /((x(i)-lambda(2,n))^2+lambda(3,n))^2;
        J(i,3)=-lambda(1,n)/((x(i)-lambda(2,n))^2+lambda(3,n))^2;
        b(i)=y(i)-F;
    end
    % Take into account error bar
    b=b./s;
    for i=1:M
       J(:,i)=J(:,i)./s';
    end
   
    % Solve the equation
    n=n+1;
    A=J'*J;
    c=J'*b';
    d=A\c;
    
    % update parameter values
    lambda(:,n)=lambda(:,n-1)+alpha*d;
    % evaluate chi sqaure
    chi2(n)=b*b';    

    % if chi square goes up stop
    if chi2(n)-chi2(n-1) > 0
        found=true;
    end
end

subplot(1,2,1)
z=[-4:0.1:6];
L=size(z,2);
for i=1:L
w(i)=lambda(1,n)/((z(i)-lambda(2,n))^2+lambda(3,n));
end
p=plot(z,w);
set(p,'linewidth',2,'color','black')
hold on
r=errorbar(x,y,s,'o');
set(r,'linewidth',2,'color','red')
xlabel('x','fontsize',14)
ylabel('f(x)','fontsize',14)
hold off

subplot(1,2,2)
r=semilogy([1:n],chi2);
set(r,'linewidth',2,'color','black')
xlabel('iteration','fontsize',14)
ylabel(texlabel('chi^2'),'fontsize',14)

