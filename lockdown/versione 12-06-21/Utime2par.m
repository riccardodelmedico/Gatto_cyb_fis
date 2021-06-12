function Utime = Utime2par(par,t)

L1 = par(1);
k1 = par(2);
t1 = par(3);
L2 = par(4);
k2 = par(5);
t2 = par(6);

Utime = L1./(1+exp(-k1.*(t'-t1))) - L2./(1+exp(-k2.*(t'-t2))); 