%**************************************************************************
%*     Section 16.4.1                                                     *
%*     filename: ch16pr04.m                                               *
%*     program listing number: 16.4                                       *
%*                                                                        *
%*     This program simulates the two-dimensional diffusion limited       *
%*     aggregates.                                                        *
%*                                                                        *
%*     Programed by Ryoichi Kawai for Computational Physics Course.       *
%*     Last modification:  03/04/2017.                                    *
%**************************************************************************
clear all
close all

% number of particles
N=20000;

L=601;  % size of the square space
L0=301; % center of the square

% set array size
x=zeros(N,1);
y=zeros(N,1);
A=zeros(L,L);

R=5;         % inner circle
R_max=3*R;   % outer circle

% seed particle
A(L0,L0)=1;
x(1)=L0;
y(1)=L0;
viscircles([x(1),y(1)],0.5,'Color','r');
axis equal
hold on

for n=2:N
    % random point on the inner circle
    theta=rand(1)*2*pi;
    x(n)=round(R*cos(theta))+L0;
    y(n)=round(R*sin(theta))+L0;
    
    % diffusion
    found=false; 
    while not(found)
        p=rand(1);
        if p<1/4
            x(n)=x(n)+1;
        elseif p<1/2
            x(n)=x(n)-1;
        elseif p<3/4
            y(n)=y(n)+1;
        else
            y(n)=y(n)-1;
        end
        r=sqrt((x(n)-L0)^2+(y(n)-L0)^2);
        if r>R_max  % out of bound - restart
            theta=rand(1)*2*pi;
            x(n)=round(R*cos(theta))+L0;
            y(n)=round(R*sin(theta))+L0;
        elseif r<R
            % hit the cluster?
            if A(x(n)+1,y(n))+A(x(n)-1,y(n))...
              +A(x(n),y(n)+1)+A(x(n),y(n)-1)>0
                found=true;
                A(x(n),y(n))=1;
                viscircles([x(n),y(n)],0.5,'Color','r');
                drawnow
                R = max(R, r+5); % adjust inner circle radius
                R_max = 3*R;     % adjust outer circle radius
                if R>L0
                    xlabel('x','fontsize',14)
                    ylabel('y','fontsize',14)
                    hold off
                    error('Out of Range')
                end
            end
        end
    end
end


     
    