%**************************************************************************
%*     Example  4.2                                                       *
%*     filename: ch04pr02.m                                               *
%*     program listing number: 4.2                                        *
%*                                                                        *
%*     This program finds roots of a cubic equation                       *
%*               a*x^3 + b*x^2 + c*x + d = 0                              *
%*                                                                        *
%*     Programed by Ryoichi Kawai for Computational Physics Course.       *
%*     Last modification:  10/13/2013.                                    *
%**************************************************************************
clear all;

% coefficients
b=-9; c=23; d=-15;

% formula for cubic polynomials
F=(3*c-b^2)/3;
G=(2*b^3 - 9*b*c + 27*d)/27;
H=G^2/4 + F^3/27;

yes = false;

if H>0
    S= (sqrt(H)-G/2)^(1./3.);
    U=-(sqrt(H)+G/2)^(1./3.);
    x = S+U-b/3;
    fprintf('Answer = %.5f\n',x);

else
    I=sqrt(G^2/4-H);
    J=I^(1./3.);
    K=acos(-G/(2*I));
    M=cos(K/3);
    N=sqrt(3)*sin(K/3);
    P=-b/3;    
    x(1) = P-J*(M+N);
    x(2) = P-J*(M-N);
    x(3) = P+2*J*M;
    fprintf('Answer = %.5f,  %.5f,  %.5f \n',x);
end    