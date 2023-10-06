%**************************************************************************
%*     Section 15.5.3                                                     *
%*     filename: ch15pr07.m                                               *
%*     program listing number: 15.7                                       *
%*                                                                        *
%*     This program simulatesa the surface growth using ballistic         *
%*     deposit model.                                                     *
%*                                                                        *
%*     Programed by Ryoichi Kawai for Computational Physics Course.       *
%*     Last modification:  02/26/2017.                                    *
%**************************************************************************
clear all  % clear all variables
close all  % close alll figures

N=20000; % number of particles to be deposited
L=200; % surface area

y=zeros(L,1); % reset the height of the surface
x=ceil(rand(N,1)*L);  % horizontal position of the particles

k=0;
for i=1:N
    y(x(i))=y(x(i))+1;  % ballistic growth
    
    % record the evolution of the growth after every 10 particles is
    % deposited
    if mod(i,10)==0
        k=k+1;
        z(k)=sum(y,1)/L;  % mean height
        w(k)=sqrt(sum((y-z(k)).^2)/L);  % roughness
    end
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

% theoretical height distribution (Gaussian formula)
mu=real(N)/L;
sg=real(N)*(L-1)/L^2;
g=zeros(n2-n1+1,1);
for n=n1:n2
    g(n-n1+1) = exp(-(n-mu)^2/(2*sg))/sqrt(2*pi*sg);
end

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
legend('simulation','theory')
legend('location','northwest')
hold off

% Figure 3: Heifht distribution
figure
bar([n1:n2],h);
hold on
p=plot([n1:n2],g);
set(p,'color','red','linewidth',2);
xlabel('height','fontsize',14)
ylabel('P(h)','fontsize',14)
legend('simulation','theory')
legend('location','northeast')