%*********************************************************************
%*     Example  1.7                                                  *
%*     filename: ch01pr01.m                                          *
%*     program listing number: 1.1                                   *
%*                                                                   *
%*     This program finds a machine epsilon by evaluating            *
%*                                                                   *
%*           1 + 2^(-n) > 1                                          *
%*                                                                   *
%*     At a certain positive n, this inenqualty becomes false.       *
%*     Then, the machine epsilon is 2^(n-1).                         *
%*                                                                   *
%*     Programed by Ryoichi Kawai for Computational Physics Course   *
%*     Revised on 10/13/2013                                         *
%*********************************************************************

clear all;

%* Find the single precision machine epsilon
epsilon = 1.0;  % create a real64 variable
n = 0;          % reset a counter

%* Reduce the value of epsilon until it becomes too small
while 1+epsilon > 1
   epsilon = epsilon/2;
   n = n+1;
end

%* The smallest single floating value which can be added to one.
epsilon = epsilon+epsilon;

%* Show the results
fprintf('\nMachine epsilon for 64 bit floating point\n');
fprintf('Stopped after %3d iterations \n',n);
fprintf('machine epsilon by computation = %16.7e \n',epsilon);
fprintf('machine epsilon by MATLAB      = %16.7e \n',eps(single(1.0)));
fprintf('1 + epsilon   = %16.8e \n',1+epsilon);
fprintf('1 + epsilon/2 = %16.8e \n',1+epsilon/2);
