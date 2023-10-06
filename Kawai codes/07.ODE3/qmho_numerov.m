function [x,psi] = qmho_numerov(E,xmax,h)

W = @(x) -(x^2-E);

N=round(xmax/h);
xmax=h*N;
x = linspace(-xmax,0.0,N+1);
w = -x.^2+E;
psi = zeros(N+1);
psi(1)=0;
psi(2)=0.001;

% shoot out to x=0  by the Numerov method
  for n=2:N
    psi(n+1) = 2*(1-5*h^2*w(n)/12)*psi(n) - (1+h^2*w(n-1)/12)*psi(n-1);
    psi(n+1) = psi(n+1)/(1+h^2*w(n+1)/12);
  end
end
