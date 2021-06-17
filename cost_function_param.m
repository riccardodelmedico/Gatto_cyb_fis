function cost_function_value = cost_function_param(Utry)

global lambda deltaE deltaP sig eta gammaI alfaI gammaA zeta gammaH alfaH ...
    gammaQ betaP betaA betaI N x0 eff1 eff2 ef1 prima_d seconda_d ...
    Lvect teta N_ott r ts xi w t_ott

%costruisco Lvect di conseguenza al controllo
Lvect = [];
for i = 1:length(Utry)
Lvect = [Lvect, Utry(i)*ones(14)];
end


% options = odeset('RelTol',1e-5,'AbsTol',1e-5);
% t = 0:1:N_ott-1;
x = ode4(@gatto_vaccini_unico,0,1,N_ott-1,x0');
tempo = 0:1:N_ott;
L = Lvect(fix(tempo)+1);

soluzione= zeros(N_ott, 24);
for j=1:1:165
    for i= 1:1:24
        soluzione(j,i)=x(24*(j-1)+i);
    end
end

S = soluzione(:,1); %suscettibili
E = soluzione(:,2); %esposti
P = soluzione(:,3); %presintomatici
I = soluzione(:,4); %infetti
A = soluzione(:,5); %asintomatici
H = soluzione(:,6); %ospedalizzati
Q = soluzione(:,7);%quarantena casalinga
D = soluzione(:,9);
S1 = soluzione(:,10); %suscettibili
E1 = soluzione(:,11); %esposti
P1 = soluzione(:,12); %presintomatici
I1 = soluzione(:,13); %infetti
A1 = soluzione(:,14); %asintomatici
H1 = soluzione(:,15); %ospedalizzati
Q1 = soluzione(:,16); %quarantena casalinga
D1 = soluzione(:,18);
S2 = soluzione(:,19); %suscettibili
E2 = soluzione(:,20); %esposti gruppo
P2 = soluzione(:,21); %presintomatici
I2 = soluzione(:,22); %infetti
A2 = soluzione(:,23); %asintomatici
arg1=S+E+P+I+A+S1+E1+P1+I1+A1+S2+E2+P2+I2+A2; %che va come costo di controllo,...
                                              %%%poichè si impedisce a questa gente di lavroare


%arg2= H+Q+H1+Q1; %questo è il  costo dovuto alla gente isolata, 
%%%ossia ocme perdita del PIl dovuto al fatto che abbiano contratto la
%%%malattia causa ''lockdwon troppo leggero'', nel senso che è il termine
%%%con cui entra in competizione il primo argomento
 
arg3 = D+D1; %come l'argomento 2, solo che in questo caso aggiugiamo il costo xi 


cost_function_value = (sum( exp(-(r).*t_ott').* ...
    (w.*L.*arg1) +...%QUESTO termine genera una matrice DA CORREGGERE
     + arg3.*(w/r + xi)) );

