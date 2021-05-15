%% Gatto model con vaccini
function xdot= gatto_vaccini(t,x)

global lambda deltaE deltaP sigm eta gammaI alfaI gammaA zeta gammaH alfaH ...
    gammaQ gammaA betaP betaA betaI N x0 eff1 eff2 ef1 prima_dose_ seconda_dose_
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
prima_dos=prima_dose_(fix(t)+1);
seconda_dos=seconda_dose_(fix(t)+1);


lambda= (betaP*P + betaI*I + betaA*A)/(S + E + P + I + A + R);
xdot = zeros(9,1);

    xdot(1) = - lambda*S-eff1*prima_dos;
    xdot(2) = lambda*S - deltaE*E;
    xdot(3) = deltaE*E - deltaP*P;
    xdot(4) = sigm*deltaP*P - (eta + gammaI + alfaI)*I;
    xdot(5) = (1 - sigm)*deltaP*P - gammaA*A;
    xdot(6) = (1 - zeta)*eta*I - (gammaH + alfaH)*H;
    xdot(7) = zeta*eta*I - gammaQ*Q;
    xdot(8) = gammaI*I + gammaA*A + gammaH*H;
    xdot(9) = alfaI*I + alfaH*H;

end