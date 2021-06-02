function cost_function_value = cost_function_param(Utry)

global lambda deltaE deltaP sig eta gammaI alfaI gammaA zeta gammaH alfaH ...
    gammaQ betaP betaA betaI N x0 eff1 eff2 ef1 prima_d seconda_d ...
    Lvect teta

%costruisco Lvect di conseguenza al controllo
Lvect = [];
for i = 1:length(Utry)
Lvect = [Lvect, Utry(i)*ones(14)];
end


options = odeset('RelTol',1e-8,'AbsTol',1e-8);
t = 0:1:N-1;
[tempo,x] = ode45('gatto_vaccini_unico',t,x0, options);

L = Lvect(fix(tempo)+1);

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
arg1=S+E+P+I+A+H+Q+S1+E1+P1+I1+A1+H1+Q1+S2+E2+P2+I2+A2;
arg2= H+Q+H1+Q1; % si suppone che esposti, asint e infetti possano effettuare il loro 
                    % lavoro in modalità remota, non andando a grave sul
                    % costo totale
cost_function_value = sum( exp(-(r).*t').* ...
    (w.*L.*(ts.*(arg1) +1 -ts) +...
    arg2.*(w/r + xi)) );
