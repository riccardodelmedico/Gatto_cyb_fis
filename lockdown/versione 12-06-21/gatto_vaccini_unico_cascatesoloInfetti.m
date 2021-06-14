function xdot= gatto_vaccini_unico_cascatesoloInfetti(t,x)

global lambda deltaE deltaP sig eta gammaI alfaI gammaA zeta gammaH alfaH ...
    gammaQ betaP betaA betaI N x0 eff1 eff2 ef1 prima_d seconda_d ...
    Lvect teta

S = x(1); %suscettibili
E = x(2); %esposti
P = x(3); %presintomatici
I1 = x(4);%infetti
I2 = x(5);
I3 = x(6);
A = x(7); %asintomatici
H = x(8); %ospedalizzati
Q = x(9); %quarantena casalinga
R = x(10); %recuperati
D = x(11); %dead
% parameters_vaccini;
% dati_vaccini;
S1 = x(12); %suscettibili
E1 = x(13); %esposti
P1 = x(14); %presintomatici
I11 = x(15);%infetti
I12 = x(16);
I13 = x(17);
A1 = x(18); %asintomatici
H1 = x(19); %ospedalizzati
Q1 = x(20); %quarantena casalinga
R1 = x(21); %recuperati
D1 = x(22); %dead
% parameters_vaccini;
% dati_vaccini;
S2 = x(23); %suscettibili
E2 = x(24); %esposti gruppo
P2 = x(25); %presintomatici
I21 = x(26);%infetti
I22 = x(27);
I23 = x(28);
A2 = x(29); %asintomatici
R2 = x(30); %recuperati

prima_dos=prima_d(fix(t)+1);
seconda_dos=seconda_d(fix(t)+1);
L=Lvect(fix(t)+1);


lambda= (betaP*(P+P1+P2) + betaI*(I3+I13+I23) + betaA*(A+A1+A2))/((S+S1+S2) + (E+E1+E2) + (I1+I2+I3+I11+I12+I13++I21+I22+I23) + (A+A1+A2) + (R+R1+R2));
xdot = zeros(30,1);

    xdot(1) = - lambda*S*(1-teta*L)^2 - eff1*prima_dos;
    xdot(2) = lambda*S*(1-teta*L)^2 - deltaE*E;
    xdot(3) = deltaE*E - deltaP*P;
    xdot(4) = sig*deltaP*P - (eta + gammaI + alfaI)*3*I1;
    xdot(5) = (eta + gammaI + alfaI)*3*I1 - (eta + gammaI + alfaI)*3*I2;
    xdot(6) = (eta + gammaI + alfaI)*3*I2 - (eta + gammaI + alfaI)*3*I3;
    xdot(7) = (1 - sig)*deltaP*P - gammaA*A;
    xdot(8) = (1 - zeta)*eta*I3 - (gammaH + alfaH)*H;
    xdot(9) = zeta*eta*I3 - gammaQ*Q;
    xdot(10) = gammaI*I3 + gammaA*A + gammaH*H;
    xdot(11) = alfaI*I3 + alfaH*H;
    %%%%POI LA SECONDA DINAMICA, CHE è ACCOPPIATA TRAMITE LAMBDA
    
    xdot(12) = - lambda*S1*(1-teta*L)^2+ eff1*prima_dos -eff2*seconda_dos;
    xdot(13) = lambda*S1*(1-teta*L)^2 - deltaE*E1;
    xdot(14) = deltaE*E1 - deltaP*P1;
    xdot(15) = sig*deltaP*P1 - ((1-ef1)*eta + gammaI +(1-ef1)*alfaI)*3*I11;
    xdot(16) = ((1-ef1)*eta + gammaI +(1-ef1)*alfaI)*3*I11 - ((1-ef1)*eta + gammaI +(1-ef1)*alfaI)*3*I12;
    xdot(17) = ((1-ef1)*eta + gammaI +(1-ef1)*alfaI)*3*I12 - ((1-ef1)*eta + gammaI +(1-ef1)*alfaI)*3*I13;
    xdot(18) = (1 - sig)*deltaP*P1 - gammaA*A1;
    xdot(19) = (1 - zeta)*(1-ef1)*eta*I13 - (gammaH +(1-ef1)*alfaH)*H1;
    xdot(20) = zeta*(1-ef1)*eta*I13 - gammaQ*Q1;
    xdot(21) = gammaI*I13 + gammaA*A1 + gammaH*H1;
    xdot(22) = (1-ef1)*alfaI*I13 + (1-ef1)*alfaH*H1;
    
    %%%%%E SIMILMENTE LA TERZA
    xdot(23) = - lambda*S2*(1-teta*L)^2+eff2*seconda_dos;
    xdot(24) = lambda*S2*(1-teta*L)^2 - deltaE*E2;
    xdot(25) = deltaE*E2 - deltaP*P2;
    xdot(26) = sig*deltaP*P2 - gammaI*3*I21;
    xdot(27) = gammaI*3*I21 - gammaI*3*I22;
    xdot(28) = gammaI*3*I22 - gammaI*3*I23;
    xdot(29) = (1 - sig)*deltaP*P2 - gammaA*A2;
    xdot(30) = gammaI*I23 + gammaA*A2;


end

