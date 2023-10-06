%**************************************************************************
%*     Section 11.4.3                                                     *
%*     filename: ch11pr03.m                                               *
%*     program listing number: 11.3                                       *
%*                                                                        *
%*     This program calculates the power spectrum of coupled oscillators. *
%*                                                                        *
%*     Uses MATLAB function fft()                                         *
%*                                                                        *
%*     Programed by Ryoichi Kawai for Computational Physics Course.       *
%*     Last modification:  02/18/2017.                                    *
%**************************************************************************
clear all

% system parameters
m=[2,4,3];
k=[2,4,4,2];
L=[1,2,1,2,8];
K=[[k(1)+k(2),-k(2),0];[-k(2),k(2)+k(3),-k(3)];[0,-k(3),k(3)+k(4)]];
b=[k(1)*L(1)-k(2)*L(2); k(2)*L(2)-k(3)*L(3);k(3)*L(3)+k(4)*(L(5)-L(4))];
M=[[1/m(1),0,0];[0,1/m(2),0];[0,0,1/m(3)]];

% initial conditions
x0=[5/3;4;16/3]; % at the equilibrium
v0=[1,0,0];

% control parameters
T=200;
N=4096;
dt=T/N;
dw=2*pi/T;
w=linspace(0.,dw*(N-1),N);
t=linspace(0.,dt*(N-1),N);

% First we calculate the trajectory of oscillators
x(1:3,1)=x0;
v(1:3,1)=v0;

% 1st step (Euler step)
x(:,2) = x(:,1)+v(:,1)*dt;

% Verlet algorithm
for i=3:N
    f = -K*x(:,i-1)+b;
    x(:,i)=2*x(:,i-1)-x(:,i-2)+(M*f)*dt^2;
    v(:,i-1)=(x(:,i)-x(:,i-2))/(2*dt);
end

% Now, we analyze the trajectories.
X=x(1,:)-x0(1); % deviation form the equilibrium
y=fft(X(1,:))/T; % Fourer transform
S = abs(y(1,:).^2); % power spectrum

subplot(1,2,1)
plot(t,x(1,:),t,x(2,:),t,x(3,:))
xlabel('$t$','fontsize',14)
ylabel('$x(t)$','fontsize',14)
axis([0 T 0 7])
h=legend('$x_1$','$x_2$','$x_3$');
set(h,'Interpreter','latex')

subplot(1,2,2)
r=plot(w,S);
axis([0 3 0 8])
xlabel('$\omega$','fontsize',14)
ylabel('$S(\omega)$','fontsize',14)
set(r,'linewidth',2)
hold on
q=plot( [2.056127, 2.056127], [0, 8], '--',...
      [1.540700, 1.540700], [0, 8], '--',...
      [0.631338, 0.631338], [0, 8], '--');
set(q,'color','black')
hold off