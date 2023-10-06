%**************************************************************************
%*     Section 15.5.3                                                     *
%*     filename: ch15pr08.m                                               *
%*     program listing number: 15.8                                       *
%*                                                                        *
%*     This program simulatesa the surface growth using ballistic         *
%*     deposit model with surface relaxation.                             *
%*                                                                        *
%*     Programed by Ryoichi Kawai for Computational Physics Course.       *
%*     Last modification:  02/26/2017.                                    *
%**************************************************************************
clear all
close all;

% To show the real time growth, set true in th following line
movie = true;

N=20000;
L=200;

y=zeros(L,1);
x=floor(rand(N,1)*L);  % random position

k=0;
for i=1:N

    j0=mod(x(i),L)+1;  % random deposition
    
    % lateral diffusion
    found = false;
    while not(found)
        j1=mod(j0-2,L)+1;  % left neighbor
        j2=mod(j0,L)+1;    % right neighbor
        
        if y(j0)<=y(j1) && y(j0)<=y(j2)  % both sides are higher
            y(j0)=y(j0)+1;               % no diffusion
            found = true;
        elseif  y(j0)> y(j1) && y(j0)>y(j2)  % both sides are lower
            if rand() > 0.5
                j0=j1;                       % diffuse to the left
            else
                j0=j2;                       % diffuse to the right
            end
        elseif y(j0)<=y(j1)                  % left side is higher
            j0=j2;                           % diffuse to the right
        else                                 % right side is higher
            j0=j1;                           % diffuse to the eft
        end
    end
    
    % record the evolution of the growth after every 10 particles is
    % deposited    if mod(i,10)==0
        k=k+1;
        z(k)=sum(y,1)/L;
        w(k)=sqrt(sum((y-z(k)).^2)/L);
end


% calculate the height distribution
n1=min(y);  % lowest
n2=max(y);  % heighest
h=zeros(n2-n1+1,1);
for i=1:L
    n=y(i)-n1+1;
    h(n)=h(n)+1;
end
h=h/sum(h);

% Figure 1: profile of the surface
b=bar([1:L],y);
hold on
plot([0,L],[0,0]);  % draw the base line
axis equal % fix the aspect ratio (needed for movie)
axis([0 L+1 0 N/L*2])  % fix the axis range
hold on
p=plot([0,L],[z(k),z(k)]);
set(p,'color','red','linewidth',1)
legend(p,'Mean height')
xlabel('x','fontsize',14)
ylabel('height','fontsize',14)
hold off

% Figure 2: Evolution of the surface roughness
figure  
p=plot(z,w);
hold on
set(p,'linewidth',2)
q=plot(z,sqrt(z));
set(q,'color','red','linewidth',2)
xlabel('height','fontsize',14)
ylabel('surface roughness','fontsize',14)
legend('simulation','ballistic')
legend('location','northwest')
hold off

% Figure 3: Heifht distribution
figure
bar([n1:n2],h);
hold on
xlabel('height','fontsize',14)
ylabel('P(h)','fontsize',14)

