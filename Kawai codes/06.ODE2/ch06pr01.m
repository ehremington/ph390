%**************************************************************************
%*     Example  6.1                                                       *
%*     filename: ch06pr01.m                                               *
%*     program listing number: 6.1-1                                      *
%*                                                                        *
%*     This program determines a launching speed that a rocket necessary  *
%*     to reach height yf in travel time tf.                              *
%*     Use function: rocket_trajectory(v,t)                               *
%*                                                                        *
%*     Programed by Ryoichi Kawai for Computational Physics Course.       *
%*     Last modification:  10/13/2013.                                    *
%**************************************************************************
clear all;

% set the boundary conditions
yf=100; tf=2;

% tolerance
tol=1e-8;

% control variable 
found = false;

% first guess
n=1;
v(n) = 50;
[y(n), vf] = rocket_trajectory(v(n),tf);
if abs(y(n)-yf) < tol
    found = true;
    v0 = v(n);
end

%second guess
n=n+1;
v(n) = 51;
[y(n), vf] = rocket_trajectory(v(n),tf);
if abs(y(n)-yf) < tol
    found = true;
    v0 = v(n);
end

% secant iteration
while not(found)
    v(n+1) = v(n) - (v(n)-v(n-1))/(y(n)-y(n-1))*(y(n)-yf);
    [y(n+1), vf] = rocket_trajectory(v(n+1),tf);
    if abs(y(n+1)-yf) < tol
       found = true;
       v0 = v(n+1);
    end
    n=n+1;
end

% show the result
fprintf('initial velocity = %.6f final velocity = %.6f \n',v(n),vf)

% plot the convergency
p=plot([1:n],v,'-o',[0,n+2],[101.9281,101.9281],'--');
xlabel('Iteration','Fontsize',14)
ylabel(texlabel('v_0'),'Fontsize',14)
axis([0 n+2 40 110])
set(p(1),'linewidth',2)
legend('numerical','exact')
legend('location','southeast')