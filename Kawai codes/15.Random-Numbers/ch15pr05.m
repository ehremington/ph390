%**************************************************************************
%*     Section 15.5.1                                                     *
%*     filename: ch15pr05.m                                               *
%*     program listing number: 15.5                                       *
%*                                                                        *
%*     This program evaluate the mean speed of the gas particles          *
%*     in a thermal equilibrium.                                          *
%*                                                                        *
%*     Programed by Ryoichi Kawai for Computational Physics Course.       *
%*     Last modification:  02/25/2017.                                    *
%**************************************************************************
clear all;

% parameters
T=300.;  % Temperature in K
k=1.380658e-23; % Boltzman constant in J/K
m=2*1.672623e-27; % H2 mass in kg

% Maxwell distribution
N=1000000;
s=sqrt(k*T/m);
v=normrnd(0.0,s,[N,3]); % 3 components (vx, vy, vz)

% speed
speed=sqrt(v(:,1).^2+v(:,2).^2+v(:,3).^2);
% mean
mean=sum(speed)/N;
% theory
exact=2*s*sqrt(2/pi);
%error
error=abs(mean-exact)/exact;

fprintf('mean speed = %10.5e (exact=%10.5e)\n',mean,exact)
fprintf('relative error = %10.5e\n',error)
