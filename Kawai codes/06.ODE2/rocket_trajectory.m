%**************************************************************************
%*     Example  6.1                                                       *
%*     filename: rocket_trajectory.m                                      *
%*     program listing number: 6.1-2                                      *
%*     Called by ch06pr01.m                                               *
%*                                                                        *
%*     This function detemines the trajectory of the rocket for a given   *
%*     initial velocity and the final time.                               *
%*        Input:  vi = initial velocity                                   *
%*                 t = final time                                         *
%*        Output:  y = final position of the rocket                       *
%*                 v = final velocity of the rocket                       *
%*                                                                        *
%*     Programed by Ryoichi Kawai for Computational Physics Course.       *
%*     Last modification:  10/13/2013.                                    *
%**************************************************************************
function [y,v]=rocket_trajectory(vi,t)

% This function calculat the height of the rocket at t;
% system parameter values
g=9.8; m=1;  C=0.01;

% control parameters
N=1000;
h=t/N;

% define force/mass as a function of v
f=@(v) -(C/m)*abs(v)*v-g;

% initial conditions
y0=0;
v0=vi;

% 4th-order Runge-Kutta
for n=1:N-1
    ky1 = v0;
    kv1 = f(v0);

    y_mid = y0 + ky1*h/2;
    v_mid = v0 + kv1*h/2;
    ky2 = v_mid;
    kv2 = f(v_mid);
      
    y_mid = y0 + ky2*h/2;
    v_mid = v0 + kv2*h/2;
    ky3 = v_mid;
    kv3 = f(v_mid);

    y_end = y0 + ky3*h;
    v_end = v0 + kv3*h;
    ky4 = v_end;
    kv4 = f(v_end);
    
    y0=y0+(ky1+2*(ky2+ky3)+ky4)*h/6;
    v0=v0+(kv1+2*(kv2+kv3)+kv4)*h/6;
end

% return the final height and velocity
y=y0;
v=v0;
end