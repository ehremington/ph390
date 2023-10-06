%**************************************************************************
%*  Example  2.2                                                          *
%*  filename: ch02pr02.m                                                  *
%*  program listing number: 2.2                                           *
%*                                                                        *
%*  This program evaluates the derivative of a given function func(x)     *
%*  at x=1 using the mean finite difference method with the accuracy      *
%*  specified by tolerance.                                               *
%*                                                                        *
%*     Programed by Ryoichi Kawai for Computational Physics Course.       *
%*     Last modification:  01/11/2015.                                    *
%**************************************************************************
clear all;

% define the function
func = @(x) x^3/3;

% Read tolerance from keyboard.
tol=input('Enter tolerance =');  % 

x=1;  % point where derivative is evaluated
h=1; % initial interval
diff_old=(func(x+h)-func(x-h))/(2*h); % initial numerical derivative

% any value bigger than tol is OK here.
delta = tol+1;

fprintf('%5s %18s %14s\n','h','derivative','error')

% Repeat until error is smaller than tolerance.
while delta>tol
    h=h/2;
    diff_new=(func(x+h)-func(x-h))/(2*h);
    delta=abs(diff_new-diff_old);
    fprintf('%10.3e %16.10e %16.10e\n',h,diff_new,delta)
    diff_old=diff_new;
end
fprintf('Tolerance is OK.\n')

