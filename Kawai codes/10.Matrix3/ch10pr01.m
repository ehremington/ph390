%**************************************************************************
%*     Example 10.1 -10.3                                                 *
%*     filename: ch10pr01.m                                               *
%*     program listing number: 10.1                                       *
%*                                                                        *
%*     This program finds eigenvalues and eigenvectos of a symmetric      *
%*     matrix using the power method.                                     *
%*                                                                        *
%*     Programed by Ryoichi Kawai for Computational Physics Course.       *
%*     Last modification:  01/08/2013.                                    *
%**************************************************************************
clear all;

% Matrix
A=[[2,1,0];[1,2,1];[0,1,2]];

% SOlution by eignvalue solver in MATLAB
[V, D] = eig(A);

% Example 10.1
[n, eig, u] = evpower(A);
fprintf('\nLargest Eigenvalue (Example 10.1)\n')
fprintf('Iteration=%d\n',n)
fprintf('Eigenvalue=%.6f (MATLAB: %.6f)\n',eig,D(3,3))
fprintf('Eigenvector=[%8.4f, %8.4f, %8.4f]\n',u);
fprintf('    (numpy):[%8.4f, %8.4f, %8.4f]\n',V(:,3));
eig_max=eig;

% Example 10.2
A=[[2,1,0];[1,2,1];[0,1,2]];
A=inv(A); 
[n, eig, u] = evpower(A);
eig=1/eig;
fprintf('\nSmallest Eigenvalue (Example 10.2)\n')
fprintf('Iteration=%d\n',n)
fprintf('Eigenvalue=%.6f (MATLAB: %.6f)\n',eig,D(1,1))
fprintf('Eigenvector=[%8.4f, %8.4f, %8.4f]\n',u);
fprintf('    (numpy):[%8.4f, %8.4f, %8.4f]\n',V(:,1));
eig_min=eig;

% Example 10.3
A=[[2,1,0];[1,2,1];[0,1,2]];
eig0=0.5*(eig_min+eig_max);
A=inv(A-eig0*eye(3,3));
[n, eig, u] = evpower(A);
eig=1.0/eig+eig0;
fprintf('\nThe Eigenvalue between the smallest and largest (Example 10.3)\n')
fprintf('Iteration=%d\n',n)
fprintf('Eigenvalue=%.6f (MATLAB: %.6f)\n',eig,D(2,2))
fprintf('Eigenvector=[%8.4f, %8.4f, %8.4f]\n',u);
fprintf('    (numpy):[%8.4f, %8.4f, %8.4f]\n',V(:,2));

