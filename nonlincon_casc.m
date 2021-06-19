function [c,ceq] = nonlincon_casc(Uing)

global lambda deltaE deltaP sig eta gammaI alfaI gammaA zeta gammaH alfaH ...
    gammaQ betaP betaA betaI N x0 eff1 eff2 ef1 prima_d seconda_d ...
    Lvect teta N_ott r ts xi w t_ott f x0_casc_inf...
    lambda deltaE deltaP sig eta gammaI alfaI gammaA zeta gammaH alfaH ...
    gammaQ betaP betaA betaI N x0 eff1 eff2 ef1 prima_d seconda_d ...
    Lvect teta
ceq = [];
Lvect = Uing;
[x_vaccini_tot33]= ode4(@gatto_vaccini_unico_cascatesoloInfetti, 0,1,N_ott-1, x0_casc_inf); 
% 
x33= zeros(N_ott, 42);
for j=1:1:N_ott
    for i= 1:1:42
        x33(j,i)=x_vaccini_tot33(42*(j-1)+i);
    end
end


c = -x33(:,1);
%così diciamo che tutte le vairabili devono essere positive
% figure(1)
% plot(c)
end