function xdot= gatto_vaccini3(t,x)

global lambda deltaE deltaP sigma eta gammaI alfaI gammaA zeta gammaH alfaH ...
    gammaQ gammaA betaP betaA betaI N x0 eff1 eff2 ef1 prima_dose_ seconda_dose_
S2 = x(1); %suscettibili
E2 = x(2); %esposti
P2 = x(3); %presintomatici
I2 = x(4); %infetti
A2 = x(5); %asintomatici
R2 = x(6); %recuperati

% parameters_vaccini;
% dati_vaccini;
%  prima_dos=prima_d(fix(t)+1);
 seconda_dos=seconda_dose_(fix(t)+1);


lambda= (betaP*P2 + betaI*I2 + betaA*A2)/(S2 + E2 + P2 + I2 + A2 + R2);
xdot = zeros(6,1);

    xdot(1) = - lambda*S2+eff2*seconda_dos;
    xdot(2) = lambda*S2 - deltaE*E2;
    xdot(3) = deltaE*E2 - deltaP*P2;
    xdot(4) = sigma*deltaP*P2 - gammaA*I2;
    xdot(5) = (1 - sigma)*deltaP*P2 - gammaA*A2;
    xdot(6) = gammaA*I2 - gammaA*A2;


end

