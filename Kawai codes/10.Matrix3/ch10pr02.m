%**************************************************************************
%*     Example 10.4                                                       *
%*     filename: ch10pr02.m                                               *
%*     program listing number: 10.2                                       *
%*                                                                        *
%*     This program finds eigenvalues and eigenvectos of a symmetric      *
%*     matrix using the inverse iteration method method.                  *
%*                                                                        *
%*     Programed by Ryoichi Kawai for Computational Physics Course.       *
%*     Last modification:  01/08/2013.                                    *
%**************************************************************************
clear all;

% Matrix
A=[[2,1,0];[1,2,1];[0,1,2]];


for i=1:6;
    % generate an initital guess
    if i<4
       q=0.5*i;
    else
       q=0.5*(i+1);
    end
    
% tolerance
tol = 1e-7;

found=false;
n=1;

% Generate a random vector
b0=rand(3,1);
b0=b0/sqrt(b0'*b0);

B=A-q*eye(3,3);

% inverse iteration method
while not(found) 
    y=gauss(B,b0);  % solve linear equation by the Gaussian elimination
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
lambda=(b1'*A*b1);
fprintf('Guess=%.3f,   Eigenvalue=%.6f \n',q,lambda)

end