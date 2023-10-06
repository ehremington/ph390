%**************************************************************************
%*     Secion 4.6.                                                        *
%*     filename: ch04pr06.m                                               *
%*     program listing number: 4.6-1                                      *
%*                                                                        *
%*     Require: rootfinding.m                                             *
%*                                                                        *
%*     This program finds trhe energy eigenvalues of a pqrticle           *
%*     in a finite square well potential using the secant root            *
%*     finding methods.                                                   *
%*                                                                        *
%*     Programed by Ryoichi Kawai for Computational Physics Course        *
%*     Revised on 01/07/2014.                                             *
%**************************************************************************

clear all;

z0=6;  %system parameter

f=@(z) z*tan(z) - sqrt(z0*z0-z*z);  %define function

%control parameters
N=100;
K=ceil(z0/pi);  %The number of roots

% initial bracket
z1=0;
z2=pi/2;

for k=1:K   
    dz = (z2-z1)/N;  % small shit
    z = rootfinding(f,z1+dz,z2-dz,10,10^(-6));
    fprintf('root=%.8f\n',z);
    z1 = z2;
    z2 = min([z0,z1 + pi]);
end
