%**************************************************************************
%*     Exercise 17.1                                                      *
%*     filename: ch17pr01.m                                               *
%*     program listing number: 17.1                                       *
%*                                                                        *
%*     This program finds the velocity distribution of one-dimensional    *
%*     ideal gas using Metropolis algorithm.                              *
%*                                                                        *
%*     Programed by Ryoichi Kawai for Computational Physics Course.       *
%*     Last modification:  03/11/2017.                                    *
%**************************************************************************
clear all
close all

N0=10000;  % number of thermalization steps
N=200000;  % number of samples
kT=1.0;    % Temperature times Boltzmann constant
m=1.0;     % mass
dv=0.1;    % maximum jump in velocity

v0=sqrt(kT/m); % thermal speed
v(1)=2.0*v0*(rand(1)-0.5); % initial velocity (radom between -0.5 and _0.5)

for i=1:N+N0-1
    found = false;
    while not(found)
        u = v(i) + dv*(2*rand(1)-1);   % candidate
        dE = m/2*(u^2-v(i)^2);         % energy change
        if exp(-dE/kT)>rand(1)            % accept or reject
            v(i+1)=u;
            found = true;
        end
    end
end

K=41;
h=histogram(v(N0+1:N+N0),K,'Normalization','pdf');
hold on
w=linspace(h.BinLimits(1),h.BinLimits(2),101);
y=1/sqrt(2*pi*kT/m) * exp(-m*w.^2/(2*kT));
% theoretical distribution (Maxwell)
p=plot(w,y);
set(p,'color','red','linewidth',2);
legend(p,'Maxwell');
hold on
legend('show')
axis([-4 4 0 0.5])
hold off
