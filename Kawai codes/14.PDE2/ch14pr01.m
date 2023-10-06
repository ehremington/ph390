%**************************************************************************
%*     Example 14.1                                                       *
%*     filename: ch14pr01.m                                               *
%*     program listing number: 14.1                                       *
%*                                                                        *
%*     This program solves 2-dimensional Laplace equation using Jacobi    *
%*     method.                                                            *
%*                                                                        *
%*     Programed by Ryoichi Kawai for Computational Physics Course.       *
%*     Last modification:  04/16/2017.                                    *
%**************************************************************************
close all;
clear all;

% parameters
a=1.0;  
b=1.0; 
V=1.0;

% spacial domain
Nx=201; % number of grids
Ny=101;
dx=0.1; % spacial step
dy=0.1;
x=linspace(-b,b,Nx);
y=linspace(0,a,Ny);

% time step
h=1./4.;  

%tolerence
tol=1.e-9;

% sampling interval
M=10;

% initial profile
phi0=zeros(Ny,Nx);  
phi0(:,1)=V;
phi0(:,Nx)=V;
phi0(1,:)=0;
phi0(Ny,:)=0;

% allocate arrays
phi1=phi0;

figure(1)
k=0;
diff=realmax();
while diff>tol
    k=k+1;
    for i=2:Ny-1
        for j=2:Nx-1
            phi1(i,j)=h*(phi0(i-1,j)+phi0(i+1,j)+phi0(i,j-1)+phi0(i,j+1));
        end
    end
    if mod(k,M)==0  % record the results
        s=(phi1-phi0).^2;
        diff=sum(s(:));
        fprintf('%d : diff=%14.6e\n',k,diff);
        pcolor(phi1); axis equal tight; shading interp;
        xlabel('x');
        ylabel('y');
        drawnow;
    end
    phi0=phi1;
end
colorbar

% Plot time evolution as cuntour
figure(2);
contour(x,y,phi1);
hold on;
[X,Y]=meshgrid(x(6:10:Nx-1),y(6:10:Ny-1));
[GX,GY]=gradient(phi1);
G=sqrt(GY.^2+GY.^2);
GX=GX./G;GY=GY./G;
quiver(X,Y,GX(6:10:Ny-1,6:10:Nx-1),GY(6:10:Ny-1,6:10:Nx-1),2)
hold off
axis equal tight;
xlabel('x');
ylabel('y');



    