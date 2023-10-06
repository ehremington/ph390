%**************************************************************************
%*     Section 7.2.1                                                      *
%*     filename: ch07pr02.m                                               *
%*     program listing number: 7.2                                        *
%*     Uses: qmho_numerov.m                                               *
%*                                                                        *
%*     This program finds an eigenvalue and eigenfunction of a quantum    *
%*     harmonic oscillator within a given bracket using the shooting      *
%*     method (Numerov and bisection methods).                            *
%*     Parity symmetry is taken into account.                             *
%*                                                                        *
%*     Programmed by Ryoichi Kawai for Computational Physics Course.      *
%*     Last modification:  10/13/2013.                                    *
%**************************************************************************
clear all;

E(1)=input('Energy Lower Blacket =');
E(2)=input('Energy Upper Blacket =');

symmetric = false;
anti_symmetric = false;
found = false;
xmax=5;
tol = 1e-8;
h = 0.001;

% Initial Lower bound
[x,psi] = qmho_numerov(E(1),xmax,h);
N=size(x,2);
error_L=(psi(N)-psi(N-1))/(x(N)-x(N-1));
error_L2=psi(N);
if abs(error_L) < tol
    found = true;
    symmetric = true;
    EM = E(1);
elseif abs(error_L2) < tol
    found = true;
    anti_symmeric = true;
    EM = E(1);
end

% Initial Upper bound
if not(found)
   [x,psi] = qmho_numerov(E(2),xmax,h);
   error_U=(psi(N)-psi(N-1))/(x(N)-x(N-1));
   error_U2=psi(N);
   if abs(error_U) < tol
      found = true;
      symmetric = true;
      EM = E(2);
   elseif abs(error_U2) < tol
      found = true;
      anti_symmeric = true;
      EM = E(2);
   end
   if error_U*error_L<0
      symmetric = true;
   elseif error_U2*error_L2<0
      anti_symmetric = true;
      error_U = error_U2;
      error_L = error_L2;
   else
      error('Blacket error!');
   end
   if symmetric & anti_symmetric
      error('Blacket error2!');
   end
end

% Begin bisection
while not(found)
    EM = sum(E)*0.5;
    [x,psi] = qmho_numerov(EM,xmax,h);
    if symmetric
        error_M=(psi(N)-psi(N-1))/(x(N)-x(N-1));
    else
        error_M=psi(N);
    end

    if abs(error_M)<tol
        found = true;
    else
        if error_M*error_L < 0
            E(2) = EM;
            error_U = error_M;
        else
            E(1) = EM;
            error_L = error_M;
        end
    end
end

% output the result
if symmetric
    fprintf('Symmetric state: ')
elseif anti_symmetric
    fprintf('Anti-Symmetric state: ')
end
fprintf('Eigenvalue= %.6f\n', EM);

X(1:N) = x(1:N);
X(N+1:2*N-1)=-x(N-1:-1:1);
if symmetric
    Y(1:N) = psi(1:N);
    Y(N+1:2*N-1)=psi(N-1:-1:1);
else
    Y(1:N) = psi(1:N);
    Y(N+1:2*N-1)=-psi(N-1:-1:1);
end

A = sum(Y(1:2:2*N-3).^2+4*Y(2:2:2*N-2).^2+Y(3:2:2*N-1).^2)*h/3;
Y = Y / sqrt(A);
p=plot(X,Y,[-xmax,xmax],[0,0],[0,0],[-1,1]);
set(p(1),'linewidth',2,'color','blue')
set(p(2:3),'color','black')
xlabel('x','fontsize',14)
ylabel(texlabel('psi(x)'),'fontsize',14)

