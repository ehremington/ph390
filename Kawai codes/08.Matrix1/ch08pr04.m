%**************************************************************************
%*     Example  8.4                                                       *
%*     filename: ch08pr04.m                                               *
%*     program listing number: 8.4                                        *
%*                                                                        *
%*     This program solves a simple linear equation with the Gaussian     *
%*     elimination with poartial pivoting and backsubstitution methods.   *
%*                                                                        *
%*     Programed by Ryoichi Kawai for Computational Physics Course.       *
%*     Last modification:  10/13/2013.                                    *
%**************************************************************************
clear all;

% Set a linear equation
N=3;
A=[[3,-1,4];[2,0,-1];[0,3,2]];
b=[2;-1;3];
P=eye(N,N);  % permutation matrix is initially an identity matrix

% Find scale factors
for i=1:N
    S(i)=max(A(i,:));
end

for n=1:N-1
    % Look for the pivot row
    j=n;
    Amax=abs(A(n,n)/S(n));
    for i=n:N
        AS=abs(A(i,n)/S(i));
        if AS > Amax
            j=i;
            Amax = AS;
        end
    end
    % Carry out pivoting
    if j ~= n
        for i=n:N
            TMP=A(n,i);
            A(n,i)=A(j,i);
            A(j,i)=TMP;
        end
        TMP=b(n);
        b(n)=b(j);
        b(j)=TMP;
        % Record the permutation
        P(n,n)=0; P(j,j)=0;
        P(n,j)=1; P(j,n)=1;
    end
    % Gaussian elimination
    for i=n+1:N
        M=-A(i,n)/A(n,n);
        A(i,n+1:N)=M*A(n,n+1:N)+A(i,n+1:N);
        b(i)=M*b(n)+b(i);
    end
    A(n+1,n)=0;  
end

% backsubstitution
for i=3:-1:1
    Ax=0;
    for j=i+1:3
        Ax = Ax+A(i,j)*x(j);
    end
    x(i) = (b(i)-Ax)/A(i,i);
end

% result
fprintf('\nA=\n')
fprintf('%8.5f  %8.5f  %8.5f\n',A')
fprintf('\nb=\n')
fprintf('%8.5f\n',b)
fprintf('\nP=\n')
fprintf('%i  %i  %i\n',P')
fprintf('\nx=\n')
fprintf('%8.5f\n',x)