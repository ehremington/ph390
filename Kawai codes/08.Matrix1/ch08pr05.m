%**************************************************************************
%*     Example  8.5                                                       *
%*     filename: ch08pr05.m                                               *
%*     program listing number: 8.5                                        *
%*                                                                        *
%*     This program calculates the inverse of a given matrix using        *
%*     Gaussi-Jordin methods.                                             *
%*                                                                        *
%*     Programed by Ryoichi Kawai for Computational Physics Course.       *
%*     Last modification:  10/13/2013.                                    *
%**************************************************************************
clear all;

% Set a linear equation
N=3;
A0=[[3,-1,4];[2,0,-1];[0,3,2]];
A=A0; % keep the original matrix
b=eye(N,N);

% scale factors
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
        TMP2(1:N) = b(n,:);
        b(n,:)=b(j,:);
        b(j,:)=TMP2(1:N);
    end
    % Gaussian elimination
    for i=n+1:N
        M=-A(i,n)/A(n,n);
        A(i,n+1:N)=M*A(n,n+1:N)+A(i,n+1:N);
        b(i,:)=M*b(n,:)+b(i,:);
    end
    A(n+1,n)=0;  
end

% backsubstitution
for i=3:-1:1
    Ax(1:N)=0;
    for j=i+1:3
        for k=1:N
            Ax(k) = Ax(k)+A(i,j)*x(j,k);
        end
    end
    for j=1:N
        x(i,j) = (b(i,j)-Ax(j))/A(i,i);
    end
end

% result
fprintf('\nInvers of A=\n')
fprintf('%8.5f  %8.5f  %8.5f\n',x)
fprintf('\nA A^(-1)=\n')
fprintf('%8.5f  %8.5f  %8.5f\n',A0*x)
