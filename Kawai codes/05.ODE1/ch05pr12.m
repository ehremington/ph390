%**************************************************************************
%*     Section 5.5.6                                                      *
%*     filename: ch05pr12.m                                               *
%*     program listing number: 5.12                                       *
%*                                                                        *
%*     This program calculate the trajectory of a particle scattered by   *
%*     Yukawa potential using he Valet method.                            *
%*                                                                        *
%*     Programed by Ryoichi Kawai for Computational Physics Course        *
%*     Revised on 01/07/2014.                                             *
%**************************************************************************
clear all;

E=input('Enter Energy =');

% control parameter
bmax=3;
N=10;
db=bmax/N;
h=0.01;

subplot(1,2,1);

for i=1:2*N+1
    
    % set impact parameter
    b(i)=db*(i-N-1);
    
    % initial conditions
    x(1)=-10; y(1)=b(i); vx(1)=sqrt(2*E); vy(1)=0;
    
    % first Euler step
    r = sqrt(x(1)^2+y(1)^2);
    Fx=x(1)/r^2 * (1/r+1) * exp(-r);
    Fy=y(1)/r^2 * (1/r+1) * exp(-r);
    x(2)=x(1)+vx(1)*h+Fx*h^2/2;
    y(2)=y(1)+vy(1)*h+Fy*h^2/2;
    
    % Verlet method
    n=2;
    while abs(x(n))<10
        r = sqrt(x(n)^2+y(n)^2);
        Fx=x(n)/r^2 * (1/r+1) * exp(-r);
        Fy=y(n)/r^2 * (1/r+1) * exp(-r);
        x(n+1)=2*x(n)-x(n-1)+Fx*h^2;
        y(n+1)=2*y(n)-y(n-1)+Fy*h^2;
        n=n+1;
    end
    
    % final velocity
    vfx=(x(n)-x(n-2))/(2*h);
    vfy=(y(n)-y(n-2))/(2*h);
    
    % scattering angle
    theta(i) = acos((vx(1)*vfx+vy(1)*vfy)/(2*E));
    
    % plot the trajctory
    plot(x(1:n),y(1:n));
    hold on
end

% draw a target atom
p=plot(0,0,'o');
set(p,'Color','red','Linewidth',3);
hold off

axis([-10,10,-10,10]);
xlabel('x');
ylabel('y');

% plot scattering angle
subplot(1,2,2);
q=plot(b,theta);
xlabel('Impact Parameter');
ylabel('Scattering Angle theta');
set(q,'Linewidth',2);


 