function cost_function_value = cost_function(Utry)

global lambda deltaE deltaP sig eta gammaI alfaI gammaA zeta gammaH alfaH ...
    gammaQ betaP betaA betaI N x0 eff1 eff2 ef1 prima_d seconda_d ...
    Lvect teta

Lvect = Utry;
options = odeset('RelTol',1e-8,'AbsTol',1e-10);
t = 0:1:N-1;
[~,x] = ode45('gatto_vaccini_unico',t,x0, options);

r=0.05;
ts=1;
xi=0;
w=65000;

S = x(:,1); %suscettibili
E = x(:,2); %esposti
P = x(:,3); %presintomatici
I = x(:,4); %infetti
A = x(:,5); %asintomatici
H = x(:,6); %ospedalizzati
Q = x(:,7); %quarantena casalinga
S1 = x(:,10); %suscettibili
E1 = x(:,11); %esposti
P1 = x(:,12); %presintomatici
I1 = x(:,13); %infetti
A1 = x(:,14); %asintomatici
H1 = x(:,15); %ospedalizzati
Q1 = x(:,16); %quarantena casalinga
S2 = x(:,19); %suscettibili
E2 = x(:,20); %esposti gruppo
P2 = x(:,21); %presintomatici
I2 = x(:,22); %infetti
A2 = x(:,23); %asintomatici

arg1=S+E+P+I+A+S1+E1+P1+I1+A1+S2+E2+P2+I2+A2; %che va come costo di controllo,
                                              %%%poichè si impedisce a questa gente di lavroare


arg2= H+Q+H1+Q1; %questo è il  costo dovuto alla gente isolata, 
%%%ossia ocme perdita del PIl dovuto al fatto che abbiano contratto la
%%%malattia causa ''lockdwon troppo leggero'', nel senso che è il termine
%%%con cui entra in competizione il primo argomento
 
arg3 = D+D1; %come l'argomento 2, solo che in questo caso aggiugiamo il costo xi 


cost_function_value = sum( exp(-(r).*t').* ...
    (w.*Lvect.*(ts.*(arg1) +1 -ts) +...
    arg2.*w/r + arg3.*(w/r + xi) ) );