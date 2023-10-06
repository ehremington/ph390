%**************************************************************************
%*     Section 5.5.7                                                      *
%*     filename: ch05pr13.m                                               *
%*     program listing number: 5.13                                       *
%*                                                                        *
%*     This program calculate the trajectory of a double pendulum by      *
%*     the 4th order Runge-Kutta method.                                  *
%*                                                                        *
%*     Programed by Ryoichi Kawai for Computational Physics Course        *
%*     Revised on 01/07/2014.                                             *
%**************************************************************************
clear all;



% system parameters 
global m1 m2 L1 L2 g
m1=2; m2=1; L1=1; L2=2; g=9.8;

% initial conditions
q1(1)=1.5;
q2(1)=3.0;
w1(1)=0;
w2(1)=0.0;
t(1)=0;
V = -(m1+m2)*g*L1*cos(q1(1))-m2*g*L2*cos(q2(1));
T = m1*L1^2*w1(1)^2/2+m2*(L1^2*w1(1)^2+L2^2*w2(1)^2 ...
    +2*L1*L2*w1(1)*w2(1)*cos(q1(1)-q2(1)))/2;
E = T+V;
 
% control parameter
tmax=200;
N=10000;
h=tmax/N;

for n=1:N-1
    % 4th-order Runge-Kutta
    dotw=DP(q1(n),q2(n),w1(n),w2(n));
    kw11=dotw(1);
    kw21=dotw(2);
    kq11=w1(n);
    kq21=w2(n);
    
    w1m = w1(n)+kw11*h/2;
    w2m = w2(n)+kw21*h/2;
    q1m = q1(n)+kq11*h/2;
    q2m = q2(n)+kq21*h/2;
    
    dotw=DP(q1m,q2m,w1m,w2m);
    kw12 = dotw(1);
    kw22 = dotw(2);
    kq12 = w1m;
    kq22 = w2m;
    
    w1m = w1(n)+kw12*h/2;
    w2m = w2(n)+kw22*h/2;
    q1m = q1(n)+kq12*h/2;
    q2m = q2(n)+kq22*h/2;

    dotw=DP(q1m,q2m,w1m,w2m);
    kw13 = dotw(1);
    kw23 = dotw(2);
    kq13 = w1m;
    kq23 = w2m;
    
    w1f = w1(n)+kw13*h;
    w2f = w2(n)+kw23*h;
    q1f = q1(n)+kq13*h;
    q2f = q2(n)+kq23*h;
    
    dotw=DP(q1f,q2f,w1f,w2f);
    kw14 = dotw(1);
    kw24 = dotw(2);
    kq14 = w1f;
    kq24 = w2f;
    
    q1(n+1)=q1(n)+(kq11+2*(kq12+kq13)+kq14)*h/6;
    q2(n+1)=q2(n)+(kq21+2*(kq22+kq23)+kq24)*h/6;
    w1(n+1)=w1(n)+(kw11+2*(kw12+kw13)+kw14)*h/6;
    w2(n+1)=w2(n)+(kw21+2*(kw22+kw23)+kw24)*h/6;

    V = -(m1+m2)*g*L1*cos(q1(n+1))-m2*g*L2*cos(q2(n+1));
    T = m1*L1^2*w1(n+1)^2/2+m2*(L1^2*w1(n+1)^2+L2^2*w2(n+1)^2 ...
        +2*L1*L2*w1(n+1)*w2(n+1)*cos(q1(n+1)-q2(n+1)))/2;
    E(n+1)=T+V;
    t(n+1)=t(1)+n*h;
end

% plot angular coordinates
subplot(1,2,1);
p=plot(t,q1,t,q2);
xlabel('t');
ylabel(texlabel('theta'));
set(p(1),'Color','blue','Linewidth',2);
set(p(2),'Color','red','Linewidth',2);
legend(p,texlabel('theta_1'),texlabel('theta_2'));
legend('Location','SouthWest');

% trajectory of the second bob in xy coordiates
subplot(1,2,2);
axis square;
x2=L1*sin(q1)+L2*sin(q2);
y2=-L1*cos(q1)-L2*cos(q2);
plot(x2,y2);
xlabel('x');
ylabel('y');