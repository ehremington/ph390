%**************************************************************************
%*     Example  2.1                                                       *
%*     filename: ch02pr01.m                                               *
%*     program listing number: 2.1                                        *
%*                                                                        *
%*  This program evaluates the derivative of a given function func(x)     *
%*  at x=1 using the three finite difference methods.                     *
%*  Errors in forward, backward and mean value methods are plotted.       *
%*                                                                        *
%*     Programed by Ryoichi Kawai for Computational Physics Course.       *
%*     Last modification:  12/25/2013.                                    *
%**************************************************************************
clear all;

% define the function
func = @(x) x^3/3;

x=1;  % the point at which derivative is evaluated.
imax=20; % number of different displacements

% title and headers
display('Absolute Errors')
fprintf('%4s %16s %18s %20s\n','h','forward','backward','mean value')

for i=1:imax
    
    % Small displacement
    h(i)=10^(-i+1);
    
    % Evaluation of numerical derivative
    d_f(i)=(func(x+h(i))-func(x))/h(i); % Forward diffrence
    d_b(i)=(func(x)-func(x-h(i)))/h(i); % Backward diffrence
    d_m(i)=(func(x+h(i))-func(x-h(i)))/(2*h(i)); % Mean value
    
    % Errors
    err_f(i)=abs(1-d_f(i));
    err_b(i)=abs(1-d_b(i));
    err_m(i)=abs(1-d_m(i));
    
    % Display the errors
    fprintf('%6.1e %18.10e %18.10e %18.10e\n',...
            h(i),err_f(i),err_b(i),err_m(i));
end

hh = h.*h; % h^2 (the error of the mean value formula is order of h^2)

% Plot data
subplot(1,2,1)  % left panel
q=semilogx(h(1:imax),d_b(1:imax),'-o',...
    h(1:imax),d_f(1:imax),'-d',...
    h(1:imax),d_m(1:imax),'-s');
title('Numerical Derivative');
xlabel('h','fontsize',14);
ylabel('Derivative','fontsize',14);
set(q(1),'Color','blue');
set(q(2),'Color','green');
set(q(3),'Color','red');    
legend(q,{'forward','backward','mean'});
legend(q,'Location','NorthWest');

subplot(1,2,2) % right panel
p=loglog(h(1:imax),err_b(1:imax),'o',...
    h(1:imax),err_f(1:imax),'d',...
    h(1:imax),err_m(1:imax),'s',...
    h(1:imax/2),h(1:imax/2),'--',h(1:imax/2),hh(1:imax/2),'--');
title('Errors in the finite difference methods');
xlabel('h','fontsize',14);
ylabel('absolute error','fontsize',14);
set(p(1),'Color','blue');
set(p(2),'Color','green');
set(p(3),'Color','red');
set(p(4),'Color','blue','LineWidth',2);
set(p(5),'Color','red','LineWidth',2);
legend(p,{'forward','backward','mean','h','h^2'});
legend(p,'Location','SouthWest');
hold off