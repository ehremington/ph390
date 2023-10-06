function [n, eig, u1]=evpower(A)

% tolerance
tol = 1e-7;
found=false;

% initial guess
x=rand(3,1);

% normalization
u0=x/sqrt((x'*x));

% power method iteration
n=1;
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

% eigenvalue
eig=u1'*A*u1;