%**************************************************************************
%*     Example 10.6                                                        *
%*     filename: ch10pr04.m                                               *
%*     program listing number: 10.4                                        *
%*                                                                        *
%*     This program finds eigenvalues and eigenvectos of a symmetric      *
%*     matrix using the QR algorithm method.                              *
%*                                                                        *
%*     Uses MATLAB function qr()                                          *
%*                                                                        *
%*     Programed by Ryoichi Kawai for Computational Physics Course.       *
%*     Last modification:  02/08/2015.                                    *
%**************************************************************************
clear all

% Define the matrix
A=[[1,-4,2];[-4,1,-2];[2,-2,-2]];

% Tolerance
tol = 1e-8;

% Magnitude of total off-diagonal elemens
S=0;
for i=1:3
    for j=i+1:3
        S=S+A(i,j)^2;
    end
end

P=eye(3,3);
n=0;
while S > tol
    n=n+1;
    [Q, R] =qr(A);  % QR decomposition
    A = R*Q; % Orthogonal transformation
    P=P*Q;   % Accumulating the transformation
    % Error evaluation
    S=0;
    for i=1:3
        for j=i+1:3
            S=S+A(i,j)^2;
        end
    end
end

fprintf('# of iterations = %d\n',n)

fprintf('\nTransformed Matrix\n')
fprintf('%8.4f  %8.4f  %8.4f\n',A')

fprintf('\nTransformation Matrix\n')
fprintf('%8.4f  %8.4f  %8.4f\n',P')           
