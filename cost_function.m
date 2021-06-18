function cost_function_value = cost_function(Utry)

global lambda deltaE deltaP sig eta gammaI alfaI gammaA zeta gammaH alfaH ...
    gammaQ betaP betaA betaI N x0 eff1 eff2 ef1 prima_d seconda_d ...
    Lvect teta N_ott r ts xi w t_ott f

%ora l'orizzonte temporale di soluzione � solo quello sui vaccini
Lvect = Utry;
tic
[x_vaccini_tot]= ode4(@gatto_vaccini_unico, 0,1,165-1, [0.99;0.01;zeros(22,1)]); 
toc

x= zeros(165, 24);
for j=1:1:165
    for i= 1:1:24
        x(j,i)=x_vaccini_tot(24*(j-1)+i);
    end
end

arg1 = x(:,1)+ x(:,2) + x(:,3) + x(:,4) + x(:,5) + x(:,10)+x(:,11)+x(:,12)...
    + x(:,13)+x(:,14)+x(:,19)+x(:,20)+ x(:,21)+x(:,22)+ x(:,23) ; %perdite in lavoro, poich� viene impedito a queste persone di lavorare

arg3 = x(:,9) + x(:,18); %costo in vite umani, sono i morti, maggiorati di un termine xi

cost_function_value = sum( exp(-(r).*t_ott').* ...
    (1-f).*(w.*Lvect.*(arg1 +1) +...
    + f.*arg3.*(w/r + xi) ) );