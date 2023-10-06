%**************************************************************************
%*     Exercise 16.3                                                      *
%*     filename: ch16pr03.m                                               *
%*     program listing number: 16.3                                       *
%*                                                                        *
%*     This program simulates the two-dimensional random walk.            *
%*                                                                        *
%*     Programed by Ryoichi Kawai for Computational Physics Course.       *
%*     Last modification:  03/04/2017.                                    *
%**************************************************************************
clear all
close all

M=100000;  % number of Brownian particles

N=1000;    % number of steps

x=zeros(M,1);  % initial position
y=zeros(M,1);
u=zeros(N,2);
w=zeros(N,2);
s=zeros(N);

for i=2:N
    g=randi([1,4],[M,1]);  % pick one of four directions to jump
    
    x(g==1)=x(g==1)+1;
    x(g==2)=x(g==2)-1;
    y(g==3)=y(g==3)+1;
    y(g==4)=y(g==4)-1;
        
    % record two sample trajectories
    u(i,1)=x(1);
    w(i,1)=y(1);
    u(i,2)=x(2);
    w(i,2)=y(2);

    % mean square displacenent
    s(i)=sum(x.^2+y.^2)/M;
end

plot(u(:,1),w(:,1),u(:,2),w(:,2));
hold on
R=sqrt(s(N));
viscircles([0,0],R,'Color','r');
axis equal

L=R*1.5;
axis([-L L -L L])
xlabel('x','fontsize',14)
ylabel('y','fontsize',14)
hold off

figure
plot(x(1:2000),y(1:2000),'.');
axis equal; % fix the aspect ratio (needed for movie)
axis([-100 100 -100 100]);  % fiz the axis range
hold on
viscircles([0,0],R,'Color','r');
xlabel('x','fontsize',14)
ylabel('y','fontsize',14)
hold off
