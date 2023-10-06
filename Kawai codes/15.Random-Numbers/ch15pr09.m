%**************************************************************************
%*     Section 15.5.3                                                     *
%*     filename: ch15pr09.m                                               *
%*     program listing number: 15.9                                       *
%*                                                                        *
%*     This program simulatesa the surface growth using ballistic         *
%*     deposit model with overhangs.                                      *
%*                                                                        *
%*     Programed by Ryoichi Kawai for Computational Physics Course.       *
%*     Last modification:  02/15/2015.                                    *
%**************************************************************************
clear all
close all

N=20000;
L=200;

% preparation for plotting
plot([0,L],[0,0]);
axis equal
axis([0 L+1 0 N/L*1.5])
hold on

y=zeros(L,1);
x=floor(rand(N,1)*L);   % random position

k=0;
for i=1:N
    j0=mod(x(i),L)+1;
    j1=mod(x(i)-1,L)+1;
    j2=mod(x(i)+1,L)+1;
    if y(j0)<y(j1) || y(j0)<y(j2)
        y(j0)=max(y(j1),y(j2));   % stick to the next site
    else
        y(j0)=y(j0)+1;            % regular deposition
    end
    
    % draw the particle
    pos = [j0-0.5 y(j0)-0.5 1. 1.];
    rectangle('Position',pos,'Curvature',[1 1],'FaceColor','Blue')
    drawnow
    
    % record the evolution of the growth after every 10 particles is
    % deposited.
    if mod(i,10)==0
        k=k+1;
        z(k)=sum(y,1)/L;
        w(k)=sqrt(sum((y-z(k)).^2)/L);
    end
end

p=plot([0,L],[z(k),z(k)]);
set(p,'color','red','linewidth',2)
legend(p,'Mean height')
hold off


