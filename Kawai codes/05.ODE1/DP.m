function [Y] = DP(q1,q2,w1,w2)

global m1 m2 L1 L2 g

A = (m1+m2)*L1;
B = m2*L2*cos(q1-q2);
C = m2*L1*cos(q1-q2);
D = m2*L2;
E = -m2*L2*w2^2*sin(q1-q2)-(m1+m2)*g*sin(q1);
F = m2*L1*w1^2*sin(q1-q2)-m2*g*sin(q2);
Y(1)=(D*E-B*F)/(A*D-B*C);
Y(2)=(A*F-C*E)/(A*D-B*C);

end


