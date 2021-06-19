function cost_function_value = cost_function_waterfall(Utry)

global lambda deltaE deltaP sig eta gammaI alfaI gammaA zeta gammaH alfaH ...
    gammaQ betaP betaA betaI N eff1 eff2 ef1 prima_d seconda_d ...
    Lvect teta N_ott r ts xi w t_ott f x0_casc_inf_opt

%ora l'orizzonte temporale di soluzione è solo quello sui vaccini
Lvect = Utry;
[x_vaccini_tot]= ode4(@gatto_vaccini_unico_cascatesoloInfetti, 0,1,N_ott-1, x0_casc_inf_opt); 

x= zeros(N_ott, 42);
for j=1:1:N_ott
    for i= 1:1:42
        x(j,i)=x_vaccini_tot(42*(j-1)+i);
    end
end

arg1 = x(:,1)+x(:,2)+x(:,3)+x(:,4)+x(:,5)+x(:,6)+x(:,7)+x(:,8)+x(:,9)+x(:,10)+x(:,11)+x(:,12)+...
    +x(:,16)+x(:,17)+x(:,18)+x(:,19)+x(:,20)+ x(:,21)+x(:,22)+ x(:,23)+x(:,24)+x(:,25)+...
    +x(:,26)+x(:,31)+x(:,32)+x(:,33)+x(:,34)+x(:,35)+x(:,36)+x(:,37)+x(:,38)+x(:,39)+...
    +x(:,40)+x(:,41); %perdite in lavoro, poiché viene impedito a queste persone di lavorare

arg2 = x(:,15) + x(:,30); %costo in vite umani, sono i morti, maggiorati di un termine xi

cost_function_value = sum( exp(-r.*t_ott').* ...
    ((1-f).*(w.*Lvect.*arg1) +...
    + f.*arg2.*(w/r + xi) ));