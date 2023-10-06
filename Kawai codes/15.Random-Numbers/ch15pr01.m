%**************************************************************************
%*     Example 15.1                                                       *
%*     filename: ch15pr01.m                                               *
%*     program listing number: 15.1                                       *
%*                                                                        *
%*     This program simulate a dice using a psueo random number           *
%*     generator.    (random() in MATLAB is not used.)                    *
%*                                                                        *
%*     Programed by Ryoichi Kawai for Computational Physics Course.       *
%*     Last modification:  02/15/2015.                                    *
%**************************************************************************
clear all

% parapmeters for random number generators
a = int64(16807);
b = int64(0); 
c = int64(2147483647);

% get a seed
x=int64(input('Seed='));

% generate uniform random numbers
N=6000;
for i=1:N
    x = mod(a*x,c);
    r(i)=double(x)/double(c);
end
 
% statistics of virtual die
P(1:6)=0;
for i=1:N
    n=ceil(6*r(i));
    P(n)=P(n)+1;
end

p=bar(P);
set(p,'facecolor',[0,0.75,0.75])
hold on
q=plot([0,7],[1000,1000],'--');
set(q,'color','red')
xlabel('face of die','fontsize',14)
ylabel('number of realization','fontsize',14)
hold off
