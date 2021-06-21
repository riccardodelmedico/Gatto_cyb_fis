function [c,ceq] = nonlincon(Uing)

global lambda deltaE deltaP sig eta gammaI alfaI gammaA zeta gammaH alfaH ...
    gammaQ betaP betaA betaI N x0 eff1 eff2 ef1 prima_d seconda_d ...
    Lvect teta N_ott r ts xi w t_ott f
ceq = [];
Lvect = Uing;
[x_vaccini_tot]= ode4(@gatto_vaccini_unico, 0,1,N_ott-1, x0'); 

x= zeros(N_ott, 24);
for j=1:1:N_ott
    for i= 1:1:24
        x(j,i)=x_vaccini_tot(24*(j-1)+i);
    end
end
c = -x(:,1); %così diciamo che tutte le vairabili devono essere positive

end