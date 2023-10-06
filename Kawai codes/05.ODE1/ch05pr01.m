%**************************************************************************
%*     Example 5.1                                                        *
%*     filename: ch05pr01.m                                               *
%*     program listing number: 5.1                                        *
%*                                                                        *
%*     This program solves Newton equation for a falling object           *
%*     using Euler and predictor-corrector methods.                       *
%*        m = mass of the object                                          *
%*        g = acceleration due to gravity                                 *
%*        gamma = frictional coefficient                                  *
%*                                                                        *
%*     Programed by Ryoichi Kawai for Computational Physics Course        *
%*     Revised on 01/07/2014.                                             *
%**************************************************************************
clear all;

% system parameters
gamma=1.0;
g=9.8;
m=1.0;

% intial condition
v_ex(1)=0.0;
v_eu(1)=0.0;
v_pc(1)=0.0;
t(1)=0.0;

% integration parameters
tmax=10; % maximum time
N=1000;  % maximum steps
h=tmax/N;% time step

for i=1:N-1
    t(i+1)=t(i)+h;  % time increment
    
    % Euler method
    F_eu=-gamma*v_eu(i)/m-g; 
    v_eu(i+1)=v_eu(i)+F_eu*h;
    
    % Predictor-Corrector method
    F_pc=-gamma*v_pc(i)/m-g;
    v_pc(i+1)=v_pc(i)+F_pc*h;  % predictor
    F_pc =-gamma/m*(v_pc(i)+v_pc(i+1))/2-g;
    v_pc(i+1)=v_pc(i)+F_pc*h;  % corrector
    
    % Exact solution
    v_ex(i+1)=m*g/gamma*(exp(-gamma*t(i+1))-1);
end

subplot(1,2,1);
q=plot(t,v_eu,t,v_pc,t,v_ex);
xlabel('t','fontsize',14);
ylabel('v(t)','fontsize',14);
set(q(1),'Color','blue','Linewidth',2);
set(q(2),'Color','red','Linewidth',2);
set(q(3),'Color','black','Linewidth',2)
legend(q,{'Euler','Predictor-Corrector','Exact'});
legend(q,'Location','NorthEast');

subplot(1,2,2);
% Plot the absolute errors
p=semilogy(t,abs(v_eu-v_ex),t,abs(v_pc-v_ex));

% Plotting options
xlabel('t','fontsize',14);
ylabel('absolute error','fontsize',14);
set(p(1),'Color','blue','Linewidth',2);
set(p(2),'Color','red','Linewidth',2);
legend(p,{'Euler','Predictor-Corrector'});
legend(p,'Location','SouthWest');