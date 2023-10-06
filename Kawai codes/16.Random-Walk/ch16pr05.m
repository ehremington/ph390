%**************************************************************************
%*     Section 16.4.2                                                     *
%*     filename: ch16pr05.m                                               *
%*     program listing number: 16.5                                       *
%*                                                                        *
%*     This program simulates the growth of dendrite.                     *
%*                                                                        *
%*     Programed by Ryoichi Kawai for Computational Physics Course.       *
%*     Last modification:  03/04/2017.                                    *
%**************************************************************************
clear all
close all

% size of the systrem
N=5000;
Lx=100;
Ly=100;

% initial setting
x=zeros(N,1);
y=zeros(N,1);
A=zeros(Lx,2*Ly);
A(:,1)=1;
H=5;
ymax=5;

% bias in y direction
e=0.01;

% initial plots
p=plot([0,Lx-1],[0.5,0.5]);
set(p,'linewidth',2,'color','black')
axis equal
axis([0 Lx 0 100])
hold on

% random position
x=floor(rand(N,1)*Lx);

% deposition process
for n=1:N
    y(n)=H;  % diffusion starts from here 
    
    found=false;
    
    % diffues in 2D space until it sticks to another.
    while not(found)
        p=rand(1);
        if p<1/4
            x(n)=mod(x(n)+1,Lx);
        elseif p<1/2
            x(n)=mod(x(n)-1,Lx);
        elseif p<3/4-e
            y(n)=y(n)+1;
            if y(n)>3*H  % if it went to high, start over.
               y(n)=H;
            end
        else
            y(n)=y(n)-1;
        end
        
        if y(n)<H
            i1=mod(x(n)-1,Lx)+1; 
            i2=mod(x(n)+1,Lx)+1;
            if A(i1,y(n)+1)+A(i2,y(n)+1)...
              +A(x(n)+1,y(n)+2)+A(x(n)+1,y(n))>0
                found=true;
                A(x(n)+1,y(n)+1)=1;
                viscircles([x(n)+1,y(n)+1],0.5,'Color','b');
                drawnow
                ymax=max(y(n),ymax); % adjust the stating height
                H=5+ymax;
                if ymax>Ly-1
                    axis([0 Lx 0 ymax*1.1])
                    xlabel('x','fontsize',14)
                    ylabel('y','fontsize',14)
                    hold off
                    error('Out of Range')
                end
            end
        end
    end
end
     
    