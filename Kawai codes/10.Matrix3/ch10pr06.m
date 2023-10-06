%**************************************************************************
%*     Section 10.5.2 - 10.5.4                                            *
%*     filename: ch10pr06.m                                               *
%*     program listing number: 10.6                                       *
%*                                                                        *
%*     This program finds energy and wavefunction of a chain of atoms     *
%*     using the QR algorithm method.                                     *
%*                                                                        *
%*     Uses MATLAB function qr()                                          *
%*                                                                        *
%*     Programed by Ryoichi Kawai for Computational Physics Course.       *
%*     Last modification:  02/09/2015.                                    *
%**************************************************************************
clear all
%
N=10;
alpha=-2;
beta=-1;

% Define the matrix
A=zeros(N,N);

for i=1:N
    A(i,i)=alpha;
end

for i=1:N-1
    j=i+1;
    A(i,j)=beta;
    A(j,i)=beta;
end   

% Tolerance
tol = 1e-8;

% Magnitude of total off-diagonal elemens
S=0;
for i=1:N
    for j=1:i-1
        S=S+A(i,j)^2;
    end
end

P=eye(N,N);
n=0;
while S > tol
    n=n+1;
    [Q, R] =qr(A);  % QR decomposition
    A = R*Q; % Orthogonal transformation
    P=P*Q;   % Accumulating the transformation
    % Error evaluation
    S=0;
    for i=1:N
        for j=1:i-1
            S=S+A(i,j)^2;
        end
    end
end

fprintf('# of iterations = %d\n',n)

fprintf('Eigenvalues\n')
for n=1:N
    E=-4*sin((N-n+1)*pi/(2*(N+1)))^2;
    fprintf('n=%2d : E=%6.4f  (%6.4f)\n',n,A(n,n),E)
end

fprintf('\nTransformation Matrix\n')
P

p=plot([1:N],P(:,1),'-o',[1:N],P(:,2),'-o',[1:N],P(:,3),'-o');
set(p(1),'linewidth',2,'color','blue')
set(p(2),'linewidth',2,'color','green')
set(p(3),'linewidth',2,'color','red')
xlabel('x','fontsize',14)
ylabel(texlabel('psi'),'fontsize',14)
legend(p,{texlabel('psi_1'),texlabel('psi_2'),texlabel('psi_3')});
legend(p,'Location','SouthEast');