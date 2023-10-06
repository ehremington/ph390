%**************************************************************************
%*     Section 19.4.2                                                     *
%*     filename: ch19pr04.m                                               *
%*     program listing number: 19.4                                       *
%*                                                                        *
%*     This program solves the Thomson problem  using genetic algorithm.  *
%*     using genetic algorithm.                                           *
%*                                                                        *
%*     Programed by Ryoichi Kawai for Computational Physics Course.       *
%*     Improved by Alex Skinner.                                          *
%*     Last modification:  04/02/2017.                                    *
%**************************************************************************
clear all;

% system parameters
N=5;
R=1;

% parameters for genetic algorithm
NP=2^12;          % population size
MP=int32(NP/2);   % half of the population
max_count=200;    % stopping condition
max_generation=10000;   % maximum generation

% parameter space (genes)
th_max=pi; th_min=0;
ph_max=2*pi; ph_min=0;

% initial parameters (genes)
th=rand(N,NP)*(th_max-th_min)+th_min;
ph=rand(N,NP)*(ph_max-ph_min)+ph_min;
th(1,:)=0;
ph(1,:)=0;
ph(2,:)=0;
th(3,1)=th_max;
ph(3,2)=ph_max;

% allocate arrays
ip=zeros(1,MP);
F=zeros(1,NP);
phs=zeros(N,NP);
ths=zeros(N,NP);

% reset the counters
found=false;
count=0;
generation=0;
Fmin=realmax();

% genetic evolution begins here
while not(found)

    generation=generation+1;
    if generation > max_generation  % too many generations
        fprintf('Max generation reached.  Terminated.\n')
        break;
    end
    
    % eveluate the initial fitness
    for i=1:NP
        F(i)=UCoulomb(th(:,i),ph(:,i),R);
    end
 
    % sort the population based on thier fitness
    [Fs, IX]=sort(F);
    ths=th(:,IX); phs=ph(:,IX); 
     
    % check if converged
    if Fs(1)>=Fmin
        count=count+1;
        if count > max_count  % exceed the waiting time limit
            found=true;
            break
        end
    else
        count=0; %reset wating time counter
        Fmin=Fs(1);
        fprintf('New low found: generation=%i, fitness= %.15f\n',generation,Fmin)
    end
    
    % find mating pairs from the survived population
    ip=randperm(MP);
    
    % generate offsprings       
    for k=1:2:NP/2
        % g ranges from -1 to 1
        g=-2*rand(2)+1;
        i=ip(k);
        j=ip(k+1);
        ths(:,k+1+MP)= ths(:,i)+g(1)*(ths(:,j)-ths(:,i));
        ths(:,k+2+MP)= ths(:,j)+g(2)*(ths(:,i)-ths(:,j));
        phs(:,k+1+MP)= phs(:,i)+g(1)*(phs(:,j)-phs(:,i));
        phs(:,k+2+MP)= phs(:,j)+g(2)*(phs(:,i)-phs(:,j));
    end
  
    % mutation rate         
    if mod(generation,10)==0
        mutation=0.8*rand();
    else
        mutation=0.3*rand();
    end
    
    % mutation  
    rm=rand(NP,1);
    for i=2:NP
        if rm(i)<mutation
            ths(:,i)=rand(N,1)*(th_max-th_min)+th_min;
            phs(:,i)=rand(N,1)*(ph_max-ph_min)+ph_min;
            ths(1,i)=0;
            phs(1:2,i)=0;
       end
    end
    
    % store new population
    th=ths; ph=phs;
end

fprintf('Lowest Energy=%.15f\n',Fmin)
for i=1:N
    fprintf('theta=%f, phi=%f\n',ths(i),phs(i))
end

for i=1:N
    X(i)=R*sin(ths(i,1))*cos(phs(i,1));
    Y(i)=R*sin(ths(i,1))*sin(phs(i,1));
    Z(i)=R*cos(ths(i,1));
end

close all;
figure(1)
plot3(X,Y,Z,'o')
axis equal;
xlim([-1,1]);
ylim([-1,1]);
zlim([-1,1]);
grid on
hold on
drawnow

bond=10.0;
for i=1:N
    for j=i+1:N
        d=sqrt((X(i)-X(j))^2+(Y(i)-Y(j))^2+(Z(i)-Z(j))^2);
        bond = min(bond,d);
    end
end
bond=bond*1.25;
for i=1:N
    for j=i+1:N
        d=sqrt((X(i)-X(j))^2+(Y(i)-Y(j))^2+(Z(i)-Z(j))^2);
        if d <= bond
            line([X(i),X(j)],[Y(i),Y(j)],[Z(i),Z(j)]);
            drawnow;
        end
    end
end
hold off