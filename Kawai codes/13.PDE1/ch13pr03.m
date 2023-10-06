%**************************************************************************
%*     Example 12.5.2                                                     *
%*     filename: ch12pr03.m                                               *
%*     program listing number: 12.3                                       *
%*                                                                        *
%*     This program solves a coupled reaction-diffusion systems based on  *
%*     the Brusselator model. The parameters are chosen to form spots.    *
%*                                                                        *
%*     Programed by Ryoichi Kawai for Computational Physics Course.       *
%*     Last modification:  02/15/2014.                                    *
%**************************************************************************
clear all
close all

% system parameters
a=2.5; b=5; % parameters for spots
Du=0.2; % diffusion constant for u
Dw=1.6; % duiffusion constant for w
L=20.0; % the size of the system (periodic boundary condition.

% control parameters
NL=100; % number of grid points
dx=L/NL; % step length
Du=Du/dx^2;
Dw=Dw/dx^2;
T=100; % total time (takes a long time to reach the final patttern)
dt=0.1/max(Du,Dw); % time step
NT=int32(T/dt);

% initial condition
u0=rand(NL,NL);
w0=rand(NL,NL);
pcolor(u0); axis equal tight; shading interp; drawnow;  

laplace_u=zeros(NL,NL);
laplace_w=zeros(NL,NL);

for k=1:NT
    t=(k-1)*dt;
    % Laplacian with periodic boundary
    laplace_u = circshift(u0,1,1)+circshift(u0,-1,1) ...
              + circshift(u0,1,2)+circshift(u0,-1,2) - 4*u0;
    laplace_w = circshift(w0,1,1)+circshift(w0,-1,1) ...
              + circshift(w0,1,2)+circshift(w0,-1,2) - 4*w0;
    % Euler step
    fu=a-(b+1)*u0+u0.^2.*w0+Du*laplace_u;
    fw=b*u0-u0.^2.*w0+Dw*laplace_w;
    u1=u0+fu*dt/2;
    w1=w0+fw*dt/2;
    
    % Laplacian at the mid time.
    laplace_u = circshift(u1,1,1)+circshift(u1,-1,1) ...
              + circshift(u1,1,2)+circshift(u1,-1,2) - 4*u1;
    laplace_w = circshift(w1,1,1)+circshift(w1,-1,1) ...
              + circshift(w1,1,2)+circshift(w1,-1,2) - 4*w1;

    % Runge-Kutta step
    fu=a-(b+1)*u1+u1.^2.*w1+Du*laplace_u;
    fw=b*u1-u1.^2.*w1+Dw*laplace_w;   
    u0=u0+fu*dt;
    w0=w0+fw*dt;
    
    pcolor(u0); axis equal tight; shading interp; drawnow;  
end
colorbar


figure(2)
pcolor(w0); axis equal tight; shading interp; drawnow;  
colorbar

