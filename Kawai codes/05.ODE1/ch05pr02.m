%**************************************************************************
%*     Example 5.2                                                        *
%*     filename: ch05pr02.m                                               *
%*     program listing number: 5.2                                        *
%*                                                                        *
%*     This program solves Newton equation for a falling object           *
%*     using Runge-Kutta 2nd and 4th order methods.                       *
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
v_rk2(1)=0.0;
v_rk4(1)=0.0;
v_pc(1)=0.0;
t(1)=0.0;

% integration parameters
tmax=10; % maximum time
N=1000;  % maximum steps
h=tmax/N;% time step

for i=1:N-1
    t(i+1)=t(i)+h;  % time increment
    
    % Runge-Kutta 2nd order
    k1=-gamma*v_rk2(i)/m-g; 
    v_mid=v_rk2(i)+k1*h/2;
    
    k2=-gamma*v_mid/m-g;
    v_rk2(i+1)=v_rk2(i)+k2*h;
    
    % RUnge-Kutta 4th order
    k1=-gamma*v_rk4(i)/m-g; 
    v_mid=v_rk4(i)+k1*h/2;
    
    k2=-gamma*v_mid/m-g;
    v_mid=v_rk4(i)+k2*h/2;

    k3=-gamma*v_mid/m-g;
    v_end=v_rk4(i)+k3*h;
    
    k4=-gamma*v_end/m-g;
    v_rk4(i+1)=v_rk4(i)+h/6*(k1+2*(k2+k3)+k4);
        
    % Exact solution
    v_ex(i+1)=m*g/gamma*(exp(-gamma*t(i+1))-1);
end

subplot(1,2,1);
q=plot(t,v_rk2,t,v_rk4,t,v_ex);
xlabel('t','fontsize',14);
ylabel('v(t)','fontsize',14);
set(q(1),'Color','blue','Linewidth',2);
set(q(2),'Color','red','Linewidth',2);
set(q(3),'Color','black','Linewidth',2)
legend(q,{'RK2','RK4','Exact'});
legend(q,'Location','NorthEast');

subplot(1,2,2);
% Plot the absolute errors
p=semilogy(t,abs(v_rk2-v_ex),t,abs(v_rk4-v_ex));

% Plotting options
xlabel('t','fontsize',14);
ylabel('absolute error','fontsize',14);
set(p(1),'Color','blue','Linewidth',2);
set(p(2),'Color','red','Linewidth',2);
legend(p,{'RK2','RK4'});
legend(p,'Location','SouthWest');