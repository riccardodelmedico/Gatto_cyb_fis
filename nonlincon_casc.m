function [c,ceq] = nonlincon_casc(Uing)

global Lvect N_ott x0_casc_inf
    
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

c = -x33(:,[1,16,31]); %le variabili dei suscettibili devono essere positive

end