%**************************************************************************
%*     Exercise 17.2                                                      *
%*     filename: ch17pr02.m                                               *
%*     program listing number: 17.2                                       *
%*                                                                        *
%*     This program simulates two-dimensional Ising model using           *
%*     the Metropolis algorithm.                                          *
%*                                                                        *
%*     Programed by Ryoichi Kawai for Computational Physics Course.       *
%*     Last modification:  03/17/2017.                                    *
%**************************************************************************
clear all
close all

rng('shuffle')

% control parameters
L=32;         % number of spins in one direction
LL=L*L;       % total number of spins
N0=20000;     % thermalization steps
N=1000000;    % number of Metropolis steps.
NS=N/1000;    % number of samples.

% temperature
T=2.0;       
beta = 1/T;

% Define arrays
kmax=int32(N/NS);
m=zeros([kmax,1]);
sigma2=zeros([kmax,1]);
E=zeros([kmax,1]);

% animation switch
movie=false;

if T>2.0
    s = int32(1-2*randi([0,1],[L,L]));   % Random (for high temperature)
else
    s = int32(ones([L,L]));              % Uniform (for low temperature)
end

% animation initial configuration
if movie
    figure
    axis([0 L+1 0 L+1])
    axis equal;
    hold on
    for i=1:L
        for j=1:L
            if s(i,j)>0
                color='blue';
            else
                color='yellow';
            end
            rectangle('Position',[i,j,1,1],'Curvature',[1 1],'FaceColor',color);
        end
    end
    drawnow;
end

% Begin Metropolis simulation 
k=0;
for n=1:N+N0
    % pick a site at random 
    i=randi([1,L]);
    j=randi([1,L]);

    % Evaluation of energy change
    i1=i+1; if i1>L;  i1=1; end;
    i2=i-1; if i2<1;  i2=L; end;
    j1=j+1; if j1>L;  j1=1; end;
    j2=j-1; if j2<1;  j2=L; end;
    s4 = s(i1,j)+s(i2,j)+s(i,j1)+s(i,j2);
    dE = 2.0*double(s4*s(i,j));

% Metropolis algorithm
    if exp(-beta*dE)> rand(1)
        s(i,j)=-s(i,j);  
        if movie
            if s(i,j)>0
                color='blue';
            else
                color='yellow';
            end
            rectangle('Position',[i,j,1,1],'Curvature',[1 1],'FaceColor',color);
            drawnow;
        end
    end
    
% measurement
    if n>N0 && mod(n,NS)==0
        k=k+1;
        m(k) = sum(s(:))/LL;                % mean magnetization
        sigma2(k) = sum(s(:).^2)/LL-m(k)^2; % variance
        % total energy
        h=0;
        for j=1:L-1
            for i=1:L-1
                 h=h+s(i,j)*(s(i+1,j)+s(i,j+1));
            end
        end
        for i=1:L
            h=h+s(i,L)*s(i,1)+s(L,i)*s(1,i);
        end
        E(k)=-h;
    end

end

% draw the final configuration
if not(movie)
    figure
    axis([0 L+1 0 L+1])
    axis equal;
    hold on
    for i=1:L
        for j=1:L
            if s(i,j)>0
                color='blue';
            else
                color='yellow';
            end
            rectangle('Position',[i,j,1,1],'Curvature',[1 1],'FaceColor',color)
        end
    end

    hold off
   drawnow;
end

% magnetization
subplot(1,2,1)
plot([1:k]*NS,m(1:k))
hold on
mu=sum(m(1:k))/k;
p=plot([NS, k*NS],[mu, mu]);
set(p,'color','red')
axis([0 N -1.1 1.1])
mx=num2str(mu,5);
legend('m(t)',['mean=' mx])
legend('location','southeast')
xlabel('steps')
ylabel('magnetization')
hold off

% energy
subplot(1,2,2)
Eavg=sum(E)/k;
C=(sum(E.^2)/k-Eavg^2)/T^2/LL;
p=plot([1:k]*NS,E(1:k)/LL);
hold on
p=plot([NS, k*NS],[Eavg/LL, Eavg/LL]);
set(p,'color','red')
axis([0 N -4 0])
mx=num2str(Eavg/LL,5);
legend('E(t)',['mean=' mx])
xlabel('steps')
ylabel('Energy/spin')
hold off

% statistics
fprintf('<m>=%.5f\n',mu)
fprintf('<E>=%.5f\n',Eavg/LL)
fprintf('<C>=%.5f\n',C)

