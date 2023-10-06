%**************************************************************************
%*     Section  8.6.3                                                     *
%*     filename: ch08pr08.m                                               *
%*     program listing number: 8.8                                        *
%*                                                                        *
%*     This program calculate the determinant of distance matrix for      *
%*     a tree graph using Gaussian elimination with partial pivoting.     *
%*                                                                        *
%*     Programed by Ryoichi Kawai for Computational Physics Course.       *
%*     Last modification:  02/01/2015.                                    *
%**************************************************************************
clear all;

% Set a linear equation
N=10;
A=[[0 , 1 , 2 , 3 , 4 , 4 , 3 , 4 , 4 , 5 ];...
   [1 , 0 , 1 , 2 , 3 , 3 , 2 , 3 , 3 , 4 ];...
   [2 , 1 , 0 , 1 , 2 , 2 , 1 , 2 , 2 , 3 ];...
   [3 , 2 , 1 , 0 , 1 , 1 , 2 , 3 , 3 , 4 ];...
   [4 , 3 , 2 , 1 , 0 , 2 , 3 , 4 , 4 , 5 ];...
   [4 , 3 , 2 , 1 , 2 , 0 , 3 , 4 , 4 , 5 ];...
   [3 , 2 , 1 , 2 , 3 , 3 , 0 , 1 , 1 , 2 ];...
   [4 , 3 , 2 , 3 , 4 , 4 , 1 , 0 , 2 , 3 ];...
   [4 , 3 , 2 , 3 , 4 , 4 , 1 , 2 , 0 , 1 ];...
   [5 , 4 , 3 , 4 , 5 , 5 , 2 , 3 , 1 , 0 ]];

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
        % Record the permutation
        P(n,n)=0; P(j,j)=0;
        P(n,j)=1; P(j,n)=1;
    end
    % Gaussian elimination
    for i=n+1:N
        M=-A(i,n)/A(n,n);
        A(i,n+1:N)=M*A(n,n+1:N)+A(i,n+1:N);
    end
    A(n+1,n)=0;  
end

p=sum(sum(P));
D=(-1)^p;
for i=1:10
    D=D*A(i,i);
end
D_GP=-(N-1)*(-2)^(N-2);
fprintf('Gaussian Elimination: %f\n',D)
fprintf('       Graham-Pollak: %f\n',D_GP)



