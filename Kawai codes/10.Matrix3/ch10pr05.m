%**************************************************************************
%*     Section 10.5.1                                                     *
%*     filename: ch10pr05.m                                               *
%*     program listing number: 10.5                                       *
%*                                                                        *
%*     This program finds eigenmodes of cpoupled harmonic oscillators.    *
%*     The higest and lowest frequencies are obgtained by the power       *
%*     method and the remaining frequency is computed by the inverse      *
%*     iteration method.                                                  *
%*                                                                        *
%*     Programed by Ryoichi Kawai for Computational Physics Course.       *
%*     Last modification:  02/08/2015.                                    *
%**************************************************************************
clear all;

% system parameters
k1=2; k2=4; k3=4; k4=2;
m1=2; m2=4; m3=3;
K=[[k1+k2, -k2, 0];[-k2,k2+k3,-k3];[0,-k3,k3+k4]];
Minv=[[1/m1,0,0];[0,1/m2,0];[0,0,1/m3]];
A0=Minv*K;

tol = 1e-7; % tolerance

% Find the largest/smallest eigenvalues
% by the power method.
for i=1:2 

    if i==1
        A=A0;  % for largest eigenvalue
    else
        A=inv(A0); % for smallest
    end

    found=false;
    n=1;
    
    x=rand(3,1); % initial guess
    u0=x/sqrt((x'*x)); % normalization

% power method iteration
    while not(found)
        x=A*x;   % update x
        u1=x/sqrt(x'*x);   % normalization
        err=sqrt( (u1-u0)'*(u1-u0) );   %error
        if err < tol
            found=true;
        else
            u0=u1;
            n=n+1;
        end
    end

    if i==1
        lambda(1)=u1'*A*u1; % largest eigenvalue
        u(1:3,1)=u1;
    else
        lambda(3)=1/(u1'*A*u1); % smallest
        u(1:3,3)=u1;
    end 
end

% The other eigenvalue by the inverse 
% iterative method.
% guess=middle between the largest and
% the lowest.
q = (lambda(1)+lambda(2))/2; 
A = A0;

b0=rand(3,1); % random vector
b0=b0/sqrt(b0'*b0);
B=A-q*eye(3,3); % shifted matrix

found=false;
n=1;

% inverse iteration method
while not(found) 
    y=gauss(B,b0);
    b1=y/sqrt(y'*y);
    if b1(1)<0 % correct the phase.
        b1=-b1;
    end
    err=sqrt( (b1-b0)'*(b1-b0) ); 
    if err < tol
        found=true;
    else
        b0=b1;
        n=n+1;
    end
end

% eigenvalue
lambda(2)=(b1'*A*b1);
u(1:3,2)=b1;

% output
for i=1:3
    fprintf('\nFrequency=%.6f \n',sqrt(lambda(i)))
    fprintf('Eigenvector\n');
    fprintf('%9.6f\n',u(:,i));
end
