function xdot= gatto_vaccini2(t,x)

global lambda deltaE deltaP sigm eta gammaI alfaI gammaA zeta gammaH alfaH ...
    gammaQ gammaA betaP betaA betaI N x0 eff1 eff2 ef1 prima_dose_ seconda_dose_
S1 = x(1); %suscettibili
E1 = x(2); %esposti
P1 = x(3); %presintomatici
I1 = x(4); %infetti
A1 = x(5); %asintomatici
H1 = x(6); %ospedalizzati
Q1 = x(7); %quarantena casalinga
R1 = x(8); %recuperati
D1 = x(9); %dead
% parameters_vaccini;
% dati_vaccini;
prima_dos=prima_dose_(fix(t)+1);
seconda_dos=seconda_dose_(fix(t)+1);


lambda= (betaP*P1 + betaI*I1 + betaA*A1)/(S1 + E1 + P1 + I1 + A1 + R1);
xdot = zeros(9,1);

    xdot(1) = - lambda*S1+eff1*prima_dos-eff2*seconda_dos;
    xdot(2) = lambda*S1 - deltaE*E1;
    xdot(3) = deltaE*E1 - deltaP*P1;
    xdot(4) = sigm*deltaP*P1 - ((1-ef1)*eta + gammaI +(1-ef1)*alfaI)*I1;
    xdot(5) = (1 - sigm)*deltaP*P1 - gammaA*A1;
    xdot(6) = (1 - zeta)*(1-ef1)*eta*I1 - (gammaH +(1-ef1)*alfaH)*H1;
    xdot(7) = zeta*(1-ef1)*eta*I1 - gammaQ*Q1;
    xdot(8) = gammaI*I1 + gammaA*A1 + gammaH*H1;
    xdot(9) = (1-ef1)*alfaI*I1 + (1-ef1)*alfaH*H1;

end