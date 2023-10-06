%**************************************************************************
%*     Example 15.2                                                       *
%*     filename: ch15pr02.m                                               *
%*     program listing number: 15.2                                       *
%*                                                                        *
%*     This program  evaluate the value of pi using the Monte Carlo       *
%*     integeration of a circle.                                          *
%*     Uses: rand() in MATLAB                                             *
%*                                                                        *
%*     Programed by Ryoichi Kawai for Computational Physics Course.       *
%*     Last modification:  02/25/2017.                                    *
%**************************************************************************
clear all
clc

% random points on a square
N=100000;
x=random('unif',-1.0,1.0,[1,N]);
y=random('unif',-1.0,1.0,[1,N]);

% points inside the circle
hit=0;
i=0;
for n=1:N
    if x(n)^2+y(n)^2 < 1
        hit=hit+1;
    end
    if mod(n,100)==0 % evaluate at every 100 
        i=i+1;
        PI(i)=hit/n*4; % estimate of pi
    end
end

p=plot([1:i]*100,PI);
hold on
q=plot([0, N], [pi,pi]);
set(q,'color','black')
xlabel('# of sampling','fontsize',14)
ylabel('V_2','fontsize',14)
axis([0 N pi*0.9 pi*1.1])
hold off

        
    