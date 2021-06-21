function [c,ceq] = nonlincon(Uing)

global x0 Lvect N_ott 
ceq = [];
Lvect = Uing;
[x_vaccini_tot]= ode4(@gatto_vaccini_unico, 0,1,N_ott-1, x0'); 

x= zeros(N_ott, 24);
for j=1:1:N_ott
    for i= 1:1:24
        x(j,i)=x_vaccini_tot(24*(j-1)+i);
    end
end
c = -x(:,[1,10,18]); %così diciamo che tutte le vairabili devono essere positive

end