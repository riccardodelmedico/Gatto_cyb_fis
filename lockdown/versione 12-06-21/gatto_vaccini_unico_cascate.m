
function xdot= gatto_vaccini_unico(t,x)

global lambda deltaE deltaP sig eta gammaI alfaI gammaA zeta gammaH alfaH ...
    gammaQ betaP betaA betaI N x0 eff1 eff2 ef1 prima_d seconda_d ...
    Lvect teta

S = x(1); %suscettibili
E1 = x(2); %esposti
E2 = x(3); %esposti2
P = x(4); %presintomatici
I = x(5); %infetti
A = x(6); %asintomatici
H1 = x(7); %ospedalizzati
H2 = x(8);
H3 = x(9);
Q = x(10); %quarantena casalinga
R = x(11); %recuperati
D = x(12); %dead
% parameters_vaccini;
% dati_vaccini;
S1 = x(13); %suscettibili
E11 = x(14); %esposti
E12 = x(15);
P1 = x(16); %presintomatici
I1 = x(17); %infetti
A1 = x(18); %asintomatici
H11 = x(19); %ospedalizzati
H12 = x(20);
H13 = x(21);
Q1 = x(22); %quarantena casalinga
R1 = x(23); %recuperati
D1 = x(24); %dead
% parameters_vaccini;
% dati_vaccini;
S2 = x(25); %suscettibili
E21 = x(26); %esposti gruppo
E22 = x(27);
P2 = x(28); %presintomatici
I2 = x(29); %infetti
A2 = x(30); %asintomatici
R2 = x(31); %recuperati

prima_dos=prima_d(fix(t)+1);
seconda_dos=seconda_d(fix(t)+1);
L=Lvect(fix(t)+1);


lambda= (betaP*(P+P1+P2) + betaI*(I+I1+I2) + betaA*(A+A1+A2))/(S + (E2+E12+E22) + (I+I1+I2) + (A+A1+A2) + R);
xdot = zeros(31,1);

    xdot(1) = - lambda*S*(1-teta*L)^2-eff1*prima_dos;
    xdot(2) = lambda*S*(1-teta*L)^2 - deltaE*2*E1;
    xdot(3) = 2*deltaE*E1 - 2*deltaE*E2; 
    xdot(4) = deltaE*E2 - deltaP*P;
    xdot(5) = sig*deltaP*P - (eta + gammaI + alfaI)*I;
    xdot(6) = (1 - sig)*deltaP*P - gammaA*A;
    xdot(7) = (1 - zeta)*eta*I - (gammaH + alfaH)*3*H1;
    xdot(8) = (gammaH + alfaH)*3*H1 - (gammaH + alfaH)*3*H2;
    xdot(9) = (gammaH + alfaH)*3*H2 - (gammaH + alfaH)*3*H3;
    xdot(10) = zeta*eta*I - gammaQ*Q;
    xdot(11) = gammaI*I + gammaA*A + gammaH*H3;
    xdot(12) = alfaI*I + alfaH*H3;
    %%%%POI LA SECONDA DINAMICA, CHE è ACCOPPIATA TRAMITE LAMBDA
    
    xdot(13) = - lambda*S1*(1-teta*L)^2+eff1*prima_dos-eff2*seconda_dos;
    xdot(14) = lambda*S1*(1-teta*L)^2 - deltaE*2*E11;
    xdot(15) = deltaE*2*E11 - deltaE*2*E12;
    xdot(16) = deltaE*E12 - deltaP*P1;
    xdot(17) = sig*deltaP*P1 - ((1-ef1)*eta + gammaI +(1-ef1)*alfaI)*I1;
    xdot(18) = (1 - sig)*deltaP*P1 - gammaA*A1;
    xdot(19) = (1 - zeta)*(1-ef1)*eta*I1 - (gammaH +(1-ef1)*alfaH)*3*H11;
    xdot(20) = (gammaH +(1-ef1)*alfaH)*3*H11 - (gammaH +(1-ef1)*alfaH)*3*H12;
    xdot(21) = (gammaH +(1-ef1)*alfaH)*3*H12 - (gammaH +(1-ef1)*alfaH)*3*H13;
    xdot(22) = zeta*(1-ef1)*eta*I1 - gammaQ*Q1;
    xdot(23) = gammaI*I1 + gammaA*A1 + gammaH*H13;
    xdot(24) = (1-ef1)*alfaI*I1 + (1-ef1)*alfaH*H13;
    
    %%%%%E SIMILMENTE LA TERZA
    xdot(25) = - lambda*S2*(1-teta*L)^2+eff2*seconda_dos;
    xdot(26) = lambda*S2*(1-teta*L)^2 - deltaE*2*E21;
    xdot(27) = deltaE*2*E21 - deltaE*2*E22;
    xdot(28) = deltaE*E22 - deltaP*P2;
    xdot(29) = sig*deltaP*P2 - gammaI*I2;
    xdot(30) = (1 - sig)*deltaP*P2 - gammaA*A2;
    xdot(31) = gammaI*I2 + gammaA*A2;


end