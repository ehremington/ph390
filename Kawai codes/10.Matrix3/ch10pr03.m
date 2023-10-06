%**************************************************************************
%*     Example 10.5                                                       *
%*     filename: ch10pr03.m                                               *
%*     program listing number: 10.3                                       *
%*                                                                        *
%*     This program finds eigenvalues and eigenvectos of a symmetric      *
%*     matrix using the Jacobi transformation method.                     *
%*                                                                        *
%*     Programed by Ryoichi Kawai for Computational Physics Course.       *
%*     Last modification:  02/08/2015.                                    *
%**************************************************************************
clear all

% Define the matrix
A=[[1,-4,2];[-4,1,-2];[2,-2,-2]];

% Tolerance
tol = 1e-4;

% Evalkuate error
S=0;
for i=1:3
    for j=i+1:3
        S=S+A(i,j)^2;
    end
end

P=eye(3,3); %Initial transformation matrix

while S > tol
    for i=1:3
        for j=i+1:3
            if A(i,j) ~= 0
                % Jacobian rotation
                if A(j,j)==A(i,i)
                    c=-1/sqrt(2);
                    s=-1/sqrt(2);
                else
                    beta=(A(j,j)-A(i,i))/(2*A(i,j));
                    t=sign(beta)/(abs(beta)+sqrt(beta^2+1));
                    c=1/sqrt(t^2+1);
                    s=t/sqrt(t^2+1);
                end
                r=s/(1+c);
                ai = c^2*A(i,i)+s^2*A(j,j)-2*s*c*A(i,j);
                aj = s^2*A(i,i)+c^2*A(j,j)+2*s*c*A(i,j);
                A(i,i)=ai;
                A(j,j)=aj;
                % Transformation matrix
                Q=eye(3,3);
                Q(i,i)=c;
                Q(j,j)=c;
                Q(i,j)=s;
                Q(j,i)=-s;
                P=P*Q;
                for k=1:3
                    if k~=i && k~=j
                        aki=c*A(k,i)-s*A(k,j);
                        akj=c*A(k,j)+s*A(k,i);
                        A(k,i)=aki;
                        A(i,k)=aki;
                        A(k,j)=akj;
                        A(j,k)=akj;
                    end
                end
                A(i,j)=0;
                A(j,i)=0;
            end
        end
    end
    
    % Evaluate error
    S=0;
    for i=1:3
        for j=i+1:3
            S=S+abs(A(i,j));
        end
    end
end

fprintf('\nTransformed Matrix\n')
fprintf('%8.4f  %8.4f  %8.4f\n',A')

fprintf('\nTransformation Matrix\n')
fprintf('%8.4f  %8.4f  %8.4f\n',P')           
