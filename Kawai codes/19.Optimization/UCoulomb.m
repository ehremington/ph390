%**************************************************************************
%*     filename: UCoulomb.m                                               *
%*                                                                        *
%*     Function called by ch19pr04.m                                      *
%**************************************************************************
function UC=UCoulomb(theta,phi,r)

N=size(theta,1);
UC=0;
X=sin(theta).*cos(phi);
Y=sin(theta).*sin(phi);
Z=cos(theta);
for i=1:N
    for j=i+1:N
        r12=sqrt((X(i)-X(j))^2+(Y(i)-Y(j))^2+(Z(i)-Z(j))^2);
        UC=UC+1.0/(r12*r);
    end
end
