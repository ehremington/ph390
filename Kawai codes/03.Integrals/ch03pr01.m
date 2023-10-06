%**************************************************************************
%*     Example  3.1                                                       *
%*     filename: ch03pr01.m                                               *
%*     program listing number: 3.1                                        *
%*                                                                        *
%*  This program integrate sin(x) from x=0 to x=pi using rectangular,     *
%*  trapezoidal and simpson methods.  Absolute errors are plotted.        *
%*                                                                        *
%*     Programed by Ryoichi Kawai for Computational Physics Course.       *
%*     Last modification:  01/06/2014.                                    *
%**************************************************************************

clear all;

% Set the lower and upper bound of the integration
a=0;
b=pi;

opt=input('Use MATLAB built-in integration? [y/n]','s');
if opt=='y'
    display('MATLAB built-in trapz will be used.');
else
    display('No MATLAB built-in function will be used.');
end

% Header of the output
display('    Absolute error in various numrical integration')
fprintf(' %3s %21s %24s %24s\n','N','Rectangular','Trapezoidal','Simpson');

% loop over different N
for k=1:10
    N=2^k;
    h(k)=(b-a)/N;
    
    % evaluate the variable and function
    for i=0:N
         x(i+1)=a+i*h(k);
         f(i+1)=sin(x(i+1));
    end
    
    % Rectangular rule
    rect=sum(f(1:N))*h(k);
    
    % Trapezoidal rule
    if opt=='y'
        % Use MATLAB trapz
        trap=trapz(x,f);
    else
        % Manual
        trap=sum(f(2:N))*h(k) + (f(1)+f(N+1))*h(k)/2;
    end
    
    % Simpson rule
    simp=(2*sum(f(1:2:N-1))+4*sum(f(2:2:N))-f(1)+f(N+1))*h(k)/3;

    % Evaluation of errors (the exat answer is 1)
    err_rect(k)=abs(1-rect);
    err_trap(k)=abs(1-trap);
    err_simp(k)=abs(1-simp);
    
    %print out the results
    fprintf('   %5d %24.16e %24.16e %24.16e \n',...
        N,err_rect(k),err_trap(k),err_simp(k));
end

% Order of the errors
h2=h.^2;
h3=h.^3;
h4=h.^4;

% Plot the results
subplot(1,1,1)
p=loglog(h,err_rect,'d',h,err_trap,'o',h,err_simp,'s',...
    h,h,'--',h,h2,'--',h,h4,'--');

% Format the plot
title('Absolute errors in the various integration methods');
xlabel('h');
ylabel('absolute error');
set(p(1),'Color','green');
set(p(2),'Color','blue');
set(p(3),'Color','red');
set(p(4),'Color','green','LineWidth',2);
set(p(5),'Color','blue','LineWidth',2);
set(p(6),'Color','red','LineWidth',2);
legend(p,{'rectangular','trapezoidal','simpson','h','h^2','h^4'});
legend(p,'Location','SouthEast');
