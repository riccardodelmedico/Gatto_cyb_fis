function xdot= gatto_vaccini_unico(t,x)

global lambda deltaE deltaP sig eta gammaI alfaI gammaA zeta gammaH alfaH ...
    gammaQ betaP betaA betaI N x0 eff1 eff2 ef1 prima_d seconda_d ...
    Lvect teta

S = x(1); %suscettibili
E = x(2); %esposti
P = x(3); %presintomatici
I = x(4); %infetti
A = x(5); %asintomatici
H = x(6); %ospedalizzati
Q = x(7); %quarantena casalinga
R = x(8); %recuperati
D = x(9); %dead
% parameters_vaccini;
% dati_vaccini;
S1 = x(10); %suscettibili
E1 = x(11); %esposti
P1 = x(12); %presintomatici
I1 = x(13); %infetti
A1 = x(14); %asintomatici
H1 = x(15); %ospedalizzati
Q1 = x(16); %quarantena casalinga
R1 = x(17); %recuperati
D1 = x(18); %dead
% parameters_vaccini;
% dati_vaccini;
S2 = x(19); %suscettibili
E2 = x(20); %esposti gruppo
P2 = x(21); %presintomatici
I2 = x(22); %infetti
A2 = x(23); %asintomatici
R2 = x(24); %recuperati

prima_dos=prima_d(fix(t)+1);
seconda_dos=seconda_d(fix(t)+1);
L=Lvect(fix(t)+1);


lambda= (betaP*(P+P1+P2) + betaI*(I+I1+I2) + betaA*(A+A1+A2))/(S + E + (I+I1+I2) + (A+A1+A2) + R);
xdot = zeros(24,1);

    xdot(1) = - lambda*S*(1-teta*L)^2-eff1*prima_dos;
    xdot(2) = lambda*S*(1-teta*L)^2 - deltaE*E;
    xdot(3) = deltaE*E - deltaP*P;
    xdot(4) = sig*deltaP*P - (eta + gammaI + alfaI)*I;
    xdot(5) = (1 - sig)*deltaP*P - gammaA*A;
    xdot(6) = (1 - zeta)*eta*I - (gammaH + alfaH)*H;
    xdot(7) = zeta*eta*I - gammaQ*Q;
    xdot(8) = gammaI*I + gammaA*A + gammaH*H;
    xdot(9) = alfaI*I + alfaH*H;
    %%%%POI LA SECONDA DINAMICA, CHE è ACCOPPIATA TRAMITE LAMBDA
    
    xdot(10) = - lambda*S1*(1-teta*L)^2+eff1*prima_dos-eff2*seconda_dos;
    xdot(11) = lambda*S1*(1-teta*L)^2 - deltaE*E1;
    xdot(12) = deltaE*E1 - deltaP*P1;
    xdot(13) = sig*deltaP*P1 - ((1-ef1)*eta + gammaI +(1-ef1)*alfaI)*I1;
    xdot(14) = (1 - sig)*deltaP*P1 - gammaA*A1;
    xdot(15) = (1 - zeta)*(1-ef1)*eta*I1 - (gammaH +(1-ef1)*alfaH)*H1;
    xdot(16) = zeta*(1-ef1)*eta*I1 - gammaQ*Q1;
    xdot(17) = gammaI*I1 + gammaA*A1 + gammaH*H1;
    xdot(18) = (1-ef1)*alfaI*I1 + (1-ef1)*alfaH*H1;
    
    %%%%%E SIMILMENTE LA TERZA
    xdot(19) = - lambda*S2*(1-teta*L)^2+eff2*seconda_dos;
    xdot(20) = lambda*S2*(1-teta*L)^2 - deltaE*E2;
    xdot(21) = deltaE*E2 - deltaP*P2;
    xdot(22) = sig*deltaP*P2 - gammaI*I2;
    xdot(23) = (1 - sig)*deltaP*P2 - gammaA*A2;
    xdot(24) = gammaI*I2 + gammaA*A2;


end

