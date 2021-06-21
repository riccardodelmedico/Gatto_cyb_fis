function xdot= gatto_vaccini_unico_cascatesoloInfetti(t,x)

global lambda deltaE deltaP sig eta gammaI alfaI gammaA zeta gammaH alfaH ...
    gammaQ betaP betaA betaI N x0 eff1 eff2 ef1 prima_d seconda_d ...
    Lvect teta

S = x(1); %suscettibili
E1 = x(2);
E2 = x(3);
E3 = x(4);%esposti
P = x(5); %presintomatici
I1 = x(6);%infetti
I2 = x(7);
I3 = x(8);
A1 = x(9);
A2 = x(10);
A3 = x(11);%asintomatici
H = x(12); %ospedalizzati
Q = x(13); %quarantena casalinga
R = x(14); %recuperati
D = x(15); %dead

S1 = x(16); %suscettibili
E11 = x(17);
E12 = x(18);
E13 = x(19);%esposti
P1 = x(20); %presintomatici
I11 = x(21);%infetti
I12 = x(22);
I13 = x(23);
A11 = x(24);
A12 = x(25);
A13 = x(26);%asintomatici
H1 = x(27); %ospedalizzati
Q1 = x(28); %quarantena casalinga
R1 = x(29); %recuperati
D1 = x(30); %dead

S2 = x(31); %suscettibili
E21 = x(32);
E22 = x(33);
E23 = x(34);%esposti gruppo
P2 = x(35); %presintomatici
I21 = x(36);%infetti
I22 = x(37);
I23 = x(38);
A21 = x(39);
A22 = x(40);
A23 = x(41);%asintomatici
R2 = x(42); %recuperati

prima_dos=prima_d(fix(t)+1);
seconda_dos=seconda_d(fix(t)+1);
L=Lvect(fix(t)+1);

lambda= (betaP*(P+P1+P2) + betaI*(I3+I13+I23) + betaA*(A3+A13+A23))/((S+S1+S2) + (E1+E2+E3+E11+E12+E13+E21+E22+E23)+...
    (I1+I2+I3+I11+I12+I13+I21+I22+I23) + (A1+A2+A3+A11+A12+A13+A21+A22+A23) + (R+R1+R2));

xdot = zeros(42,1);

    xdot(1) = - lambda*S*(1-teta*L)^2 - eff1*prima_dos;
    xdot(2) = lambda*S*(1-teta*L)^2 - deltaE*3*E1;
    xdot(3) = deltaE*3*E1 - deltaE*3*E2;
    xdot(4) = deltaE*3*E2- deltaE*3*E3;
    xdot(5) = deltaE*3*E3 - deltaP*P; %pauci
    xdot(6) = sig*deltaP*P - (eta + gammaI + alfaI)*3*I1;
    xdot(7) = (eta + gammaI + alfaI)*3*I1 - (eta + gammaI + alfaI)*3*I2;
    xdot(8) = (eta + gammaI + alfaI)*3*I2 - (eta + gammaI + alfaI)*3*I3;
    xdot(9) = (1 - sig)*deltaP*P - gammaA*3*A1;
    xdot(10) = gammaA*3*A1 - gammaA*3*A2;
    xdot(11) = gammaA*3*A2 - gammaA*3*A3;
    xdot(12) = (1 - zeta)*eta*3*I3 - (gammaH + alfaH)*H;
    xdot(13) = zeta*eta*3*I3 - gammaQ*Q;
    xdot(14) = gammaI*3*I3 + gammaA*3*A3 + gammaH*H;
    xdot(15) = alfaI*3*I3 + alfaH*H;
    %%%%POI LA SECONDA DINAMICA, CHE è ACCOPPIATA TRAMITE LAMBDA
    
    xdot(16) = - lambda*S1*(1-teta*L)^2+ eff1*prima_dos -eff2*seconda_dos;
    xdot(17) = lambda*S1*(1-teta*L)^2 - deltaE*3*E11;
    xdot(18) = deltaE*3*E11-deltaE*3*E12;
    xdot(19) = deltaE*3*E12 - deltaE*3*E13;
    xdot(20) = deltaE*3*E13 - deltaP*P1;
    xdot(21) = sig*deltaP*P1 - ((1-ef1)*eta + gammaI +(1-ef1)*alfaI)*3*I11;
    xdot(22) = ((1-ef1)*eta + gammaI +(1-ef1)*alfaI)*3*I11 - ((1-ef1)*eta + gammaI +(1-ef1)*alfaI)*3*I12;
    xdot(23) = ((1-ef1)*eta + gammaI +(1-ef1)*alfaI)*3*I12 - ((1-ef1)*eta + gammaI +(1-ef1)*alfaI)*3*I13;
    xdot(24) = (1 - sig)*deltaP*P1 - gammaA*3*A11;
    xdot(25) = gammaA*3*A11 - gammaA*3*A12;
    xdot(26) = gammaA*3*A12 - gammaA*3*A13;
    xdot(27) = (1 - zeta)*(1-ef1)*eta*3*I13 - (gammaH +(1-ef1)*alfaH)*H1;
    xdot(28) = zeta*(1-ef1)*eta*3*I13 - gammaQ*Q1;
    xdot(29) = gammaI*3*I13 + gammaA*3*A13 + gammaH*H1;
    xdot(30) = (1-ef1)*alfaI*3*I13 + (1-ef1)*alfaH*H1;
    
    %%%%%E SIMILMENTE LA TERZA
    xdot(31) = - lambda*S2*(1-teta*L)^2+eff2*seconda_dos;
    xdot(32) = lambda*S2*(1-teta*L)^2 - deltaE*3*E21;
    xdot(33) = deltaE*3*E21 - deltaE*3*E22;
    xdot(34) = deltaE*3*E22 - deltaE*3*E23;
    xdot(35) = deltaE*3*E23 - deltaP*P2;
    xdot(36) = sig*deltaP*P2 - gammaI*3*I21;
    xdot(37) = gammaI*3*I21 - gammaI*3*I22;
    xdot(38) = gammaI*3*I22 - gammaI*3*I23;
    xdot(39) = (1 - sig)*deltaP*P2 - gammaA*3*A21;
    xdot(40) = gammaA*3*A21 - gammaA*3*A22;
    xdot(41) = gammaA*3*A22 - gammaA*3*A23;
    xdot(42) = gammaI*3*I23 + gammaA*3*A23;

end

