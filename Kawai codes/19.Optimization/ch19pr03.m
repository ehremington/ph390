%**************************************************************************
%*     Section 19.4.1                                                     *
%*     filename: ch19pr03.m                                               *
%*     program listing number: 19.3                                       *
%*                                                                        *
%*     This program fits a Gaussian disitrbution to a noisy data set      *
%*     using genetic algorithm.                                           *
%*                                                                        *
%*     Programed by Ryoichi Kawai for Computational Physics Course.       *
%*     Last modification:  04/02/2017.                                    *
%**************************************************************************


% Generate an experimental data set
N=13;
x=[-1.98,-1.48,-1.00,-0.50,-0.02,0.51,1.02,1.52,1.99,2.52,3.00,3.52,3.99];
y=[ 0.10, 0.11, 0.31, 0.65, 1.28,1.79,2.13,1.92,1.73,0.70,0.31,0.14,0.14];
s=[ 0.14, 0.11, 0.15, 0.18, 0.18,0.12,0.11,0.10,0.16,0.15,0.10,0.17,0.18];

% theoretical data space
K=101;
dx=(x(N)-x(1))/(K-1);
X=linspace(x(1),x(N),K);

% control parameter for genetic algorithm
NP=2^10;
MP=int32(NP/2);
max_count=200;
max_generation=10000;

% parameter spce (genes)
a_max=5; a_min=0;
b_max=5; b_min=-5;
c_max=5; c_min=0.01;

% initial parameters (genes)
a=rand(NP,1)*(a_max-a_min)+a_min;
b=rand(NP,1)*(b_max-b_min)+b_min;
c=rand(NP,1)*(c_max-c_min)+c_min;

%  allocate arrays
ip=zeros(1,MP);
F=zeros(1,NP);
G=zeros(N,NP);

% graphics setting
close all;
movie = true;
if movie
    figure(1);
    axis([x(1) x(end) -0.5 2.5]);
end

% reset the counters
found=false;
count=0;
generation=0;
Fmin=realmax();   % some large number

% genetic evolution begins here
while not(found)
    
    generation=generation+1;
    if generation > max_generation  % too many generations
        fprintf('Max generation reached.  Terminated.\n');
        break;
    end
    
    % eveluate the initial fitness
    for i=1:NP
        for j=1:N
            G(j,i)=a(i)*exp(-(x(j)-b(i))^2/c(i));
        end
        F(i)=sum((G(:,i)-y(:)).^2./s(:))/N;
    end
 
    % sort the population based on thier fitness
    [Fs, IX]=sort(F);
    as=a(IX); bs=b(IX); cs=c(IX);

    if movie
        % plot current best fitting
        Y=as(1)*exp(-(X-bs(1)).^2/cs(1));
        r=plot(X,Y);
        set(r,'linewidth',2)
        hold on
        r=errorbar(x,y,s,'o');
        set(r,'linewidth',2,'color','red')
        xlabel('x','fontsize',14)
        ylabel('f(x)','fontsize',14)
        hold off
        axis([x(1) x(N) -0.5 2.5]);
        drawnow
        pause(0.2)
    end
    
    % check if converged
    if Fs(1)>=Fmin
        count=count+1;
        if count > max_count  % exceed the waiting time limit
            found=true;
            break
        end
    else
        count=0;            %reset wating time counter
        Fmin=Fs(1);         % new lowest fitness
        a0=as(1); b0=bs(1); c0=cs(1);
        fprintf('New low found: generation=%i, fitness= %.15f, a=%f, b=%f, c=%f\n',...
            generation,Fmin,a0,b0,c0)
    end
    
    % find mating pairs from the survived population
    ip=randperm(MP); 
    
    % generate offsprings
    for k=1:2:MP
        g=rand(6);
        i=ip(k);
        j=ip(k+1);
        as(k+1+NP/2)= as(i)+g(1)*(as(j)-as(i)+0.01);
        as(k+2+NP/2)= as(j)+g(2)*(as(i)-as(j)+0.01);
        bs(k+1+NP/2)= bs(i)+g(3)*(bs(j)-bs(i)+0.01);
        bs(k+2+NP/2)= bs(j)+g(4)*(bs(i)-bs(j)+0.01);
        cs(k+1+NP/2)= cs(i)+g(5)*(cs(j)-cs(i)+0.01);
        cs(k+2+NP/2)= cs(j)+g(6)*(cs(i)-cs(j)+0.01);
    end
 
    % mutation rate (every 10 generations, mutation burst happens
    if mod(generation,10)==0
        mutation=0.8;
    else
        mutation=0.1;
    end
    
    %  mutation
    rm=rand(1,NP);
    for i=2:NP
        if rm(i)<mutation
            as(i)=rand(1)*(a_max-a_min)+a_min;
            bs(i)=rand(1)*(b_max-b_min)+b_min;
            cs(i)=rand(1)*(c_max-c_min)+c_min;
        end
    end
    
    % store genes of new population
    a=as; b=bs; c=cs;
 
end
    
fprintf('Best fit=%.15f,  a=%f, b=%f, c=%f\n',Fmin,a0,b0,c0)

% plot the best fitting
if not(movie)
    figure(1);
    axis([x(1) x(end) -0.5 2.5])
end

Y=a0*exp(-(X-b0).^2/c0);
r=plot(X,Y);
set(r,'linewidth',2)
hold on
r=errorbar(x,y,s,'o');
set(r,'linewidth',2,'color','red')
xlabel('x','fontsize',14)
ylabel('f(x)','fontsize',14)
axis([x(1) x(N) -0.5 2.5]);
hold off
drawnow
