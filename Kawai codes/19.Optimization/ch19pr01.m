%**************************************************************************
%*     Example 19.2                                                       *
%*     filename: ch19pr01.m                                               *
%*     program listing number: 19.1                                       *
%*                                                                        *
%*     This program attempt to find a grobal minimum of a fitness         *
%*     function U(x) using simulated annealing.                                         *
%*                                                                        *
%*     Programed by Ryoichi Kawai for Computational Physics Course.       *
%*     Last modification:  04/01/2017.                                    *
%**************************************************************************
close all
clc

L=241; % size of the descretized configuration space
dx=0.1;
xmin=(121-L)*dx; % bound of the configuration space
xmax=(L-121)*dx;
x=linspace(xmin,xmax,L);

% fitness function
U=@(x) -3.*cos(2.*x)+0.2*x.^2+3.0;

% graphics parameter
movie=false;  % movie slows down the simulation

NP=2^3;  % population size
p=rand(NP,1)*(xmax-xmin)+xmin;  % initial configuration
T=5; % initial temperature
E=U(p); % initial fitness
dp=0.1; % Metropolis max step length

NT=200; % Total cooling steps
NS=2000; % Total thermalization steps
RT=0.98; % Cooling rate

if movie
    figure(1)
    set(gcf,'units','inches','position',[1,1,6,5])
end

for k=1:NT % Cooling Loop
    temp(k)=T;

    for i=1:NS % Thermalization Loop (Metropolis)
        
        for j=1:NP % Loop over population

            found=false;
            while not(found)
                p0=p(j)+(1-2*rand(1))*dp;
                E0=U(p0);
                if exp(-(E0-E(j))/T) > rand(1)
                    found=true;
                    if p0<xmin || p0>xmax
                        p(j)=rand(1)*(xmax-xmin)+xmin;
                        E(j)=U(p(j));
                    else
                        p(j)=p0;
                        E(j)=E0;
                    end
                end
            end
        end

        if movie && mod(i,10)==0 % update movie
            plot(x,U(x))
            xlim([xmin,xmax])
            ylim([-1,19])
            hold on
            for j=1:NP
                rectangle('Position',[p(j)-0.25,U(p(j))-0.25,0.5,0.5],...
                    'Curvature',[1 1],'FaceColor','b');
            end
            drawnow
            hold off
        end

    end
    Emin(k)=min(E);
    fprintf('T=%f,   Emin=%f\n',T,Emin(k))
    
    T=T*0.98; % Exponential cooling schedule
end

fprintf('Final Temperature = %d\n',temp(200));
fprintf('Final Lowest Fitness = %d\n', Emin(200));
    
figure(1)
q=plot(x,U(x));
xlim([xmin,xmax])
ylim([-1,19])
xlabel('x','fontsize',14)
ylabel('Fitness','fontsize',14)
hold on
for j=1:NP
	rectangle('Position',[p(j)-0.25,U(p(j))-0.25,0.5,0.5],'Curvature',[1 1]);
end
drawnow
hold off

figure(2)
q=plot([1:200],temp,[1:200],Emin);
set(q(1),'color','black')
set(q(2),'color','red')
ylim([-2,6]);
xlabel('Steps','fontsize',14)
ylabel('Temperature/Lowest Fitness Value','fontsize',14)
legend('Temperature','Lowest Fitness Value')
    
    