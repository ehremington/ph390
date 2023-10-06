%**************************************************************************
%*     Example 19.3                                                       *
%*     filename: ch19pr02.m                                               *
%*     program listing number: 19.2                                       *
%*                                                                        *
%*     This program findw a grobal minimum of a fitness function U(x)     *
%*     using genetic algorithm.                                           *
%*                                                                        *
%*     Programed by Ryoichi Kawai for Computational Physics Course.       *
%*     Last modification:  04/01/2017.                                    *
%**************************************************************************
close all

L=121; % size of the descretized configuration space
dx=0.1;
xmin=0; % bound of the configuration space
xmax=(L-1)*dx;
x=[0:L-1]*dx;
U=@(x) cos(5*x)-2*sin(3.5*x)+0.5*cos(x+0.5)+4; % fitness function

% parameter sof genetic algorithm
N=2^8; % size of population
p=rand(N,1)*xmax; % initial population
max_count=100; % waiting time
mutaion=0.1;  % mutation rate (fixed)

% initial fitness
f=U(p);
[f_sort,ix]=sort(f);
fmin=f_sort(1);
q=p(ix);  % q=individual in the order of their fitness



% plot initial population
figure(1)
set(gcf,'units','inches','position',[1,1,6,6])
plot(x,U(x))
axis([xmin xmax -1 11])
hold on
for i=1:N;
    rectangle('Position',[p(i)-0.2,U(p(i))-0.2,0.4,0.4],'Curvature',[1 1]);
end
drawnow
hold off

count=0;
found=false;
ip=zeros(N/2,1);
generation=0;

while not(found)
    generation=generation+1;
    pause(0.1)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % select mating pairs from survivors (Knuth shuffle)
    for i=1:N/2
        j=ceil(rand(1)*i);
        ip(i)=ip(j);
        ip(j)=i;
    end
    % The following MATLAB builtin function can be used in place for
    % the above Knuth shuffle
    % ip=randperm(N/2);  
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    % generate offspring
    for k=1:2:N/2

        i=ip(k);
        j=ip(k+1);
        g=rand(2);
       
        % replace the deads with the new borns
        q(k+1+N/2)= q(i)+g(1)*(q(j)-q(i)+0.01);  % inheritance
        q(k+2+N/2)= q(j)+g(2)*(q(i)-q(j)+0.01);  % inheritance
    end
    
    % mutation
    if mod(generation,10)==0
        mutation=0.8;
    else
        mutation=0.1;
    end
    if rand(1)<mutaion
        i=ceil(rand(1)*N);
        q(i)=rand(1)*xmax;
    end
    
    p=q; % store new population
    
    % evaluate fitness
    f=U(p);
    % plot the current generation
    plot(x,U(x))
    axis([xmin xmax -1 11])
    hold on
    for i=1:N;
        rectangle('Position',[p(i)-0.2,U(p(i))-0.2,0.4,0.4],'Curvature',[1 1]);
    end
    drawnow
    pause(0.5);
    hold off

    %  sort the population
    [f_sort,ix]=sort(f);
    q=p(ix);
    
    % check if converged
    if f_sort(1)>=fmin
        count=count+1;
        if count > max_count
            found=true;  % waited long enough
        end
    else
        count=0; % reset wating couter
        fmin=f_sort(1);
        fprintf('New low found: generation=%d fitness=%.15f\n',generation,fmin)
    end
    
    bestfit(generation)=fmin;

end

fprintf('optimal x=%f,   U(x)=%f\n',q(1),f_sort(1))

figure(2)
plot([1:generation],bestfit(1:generation),'-or')
hold on
axis([0 generation 0.4 1]) 
xlabel('Generation','fontsize',14)
ylabel('Best Finess','fontsize',14)