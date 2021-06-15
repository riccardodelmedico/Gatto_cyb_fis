function cost_function_value = cost_function(Utry)

global lambda deltaE deltaP sig eta gammaI alfaI gammaA zeta gammaH alfaH ...
    gammaQ betaP betaA betaI N x0 eff1 eff2 ef1 prima_d seconda_d ...
    Lvect teta N_ott r ts xi w t_ott

%ora l'orizzonte temporale di soluzione è solo quello sui vaccini

Lvect = Utry;
options = odeset('RelTol',1e-7,'AbsTol',1e-8);

[~,x] = ode45('gatto_vaccini_unico',t_ott,x0, options);

arg1 = x(:,1)+ x(:,2) + x(:,3) + x(:,4) + x(:,5) + x(:,10)+x(:,11)+x(:,12)...
    + x(:,13)+x(:,14)+x(:,19)+x(:,20)+ x(:,21)+x(:,22)+ x(:,23) ;

arg3 = x(:,9) + x(:,18);

%arg1=S+E+P+I+A+S1+E1+P1+I1+A1+S2+E2+P2+I2+A2; %che va come costo di controllo,
                                              %%%poichè si impedisce a questa gente di lavroare


%arg2= H+Q+H1+Q1; %questo è il  costo dovuto alla gente isolata, 
%%%ossia ocme perdita del PIl dovuto al fatto che abbiano contratto la
%%%malattia causa ''lockdwon troppo leggero'', nel senso che è il termine
%%%con cui entra in competizione il primo argomento
 
%arg3 = D+D1; %come l'argomento 2, solo che in questo caso aggiugiamo il costo xi 


cost_function_value = sum( exp(-(r).*t_ott').* ...
    (w.*Lvect.*(ts.*(arg1) +1 -ts) +...
    + arg3.*(w/r + xi) ) );