function [c,ceq] = nonlincon(Uing)

ceq = [];

c(1) = Uing(4)-Uing(1); % L2<L1 --> Uing(4)<Uing(1) --> Uing(4)-Uing(1)<0 
c(2) = Uing(3)-Uing(6); % t1<t2 --> Uing(3)<Uing(6) --> Uing(3)-Uing(6)<0
