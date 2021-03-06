global  deltaE deltaP sig eta gammaI alfaI gammaA zeta gammaH ...
    alfaH gammaQ gammaA betaP betaI betaA eff1 eff2 ef1 teta R0

deltaE= 1 / 3.32;
deltaP= 1 / 0.75;
eta= 1 / 4.05;
gammaI= 1 / 14.32;
gammaQ=gammaI;
gammaH= gammaI;
gammaA= 2*gammaI;
alfaI= 1 / 24.23;
sig= 0.25; % (1-sigma) frazione degli asintomatici
betaA_betaP= 0.033; %betaA asymptomatic transmission rate
betaI_betaA= 1.03;
% betaP1/ betaP== 0.82; % betaP1 after restriction on february 22, 2020
% betaP2/ betaP1== 0.66; % betaP2 after restriction on march 8, 2020
% deltat0= 34.94;
% omega= 7.84;
zeta= 0.4; % frazione dei sintomatici infettati in quarantena
alfaH = alfaI;
% rs= 0.5;
eff1 = 0.9; %efficacia della prima dose nel prevenire da infezione
eff2 = 1; %efficacia della seconda dose nel prevenire da infezione
ef1 = 0.7; %efficacia della prima dose nel prevenire da gravi sintomi
teta = 0.75; %efficacia del lockdown

options_lockdown = optimoptions('fmincon','Display','iter-detailed','Algorithm','active-set','FunValCheck','on','MaxFunctionEvaluation',2.7e04);

x1=1/deltaP;
x2=sig/(eta+alfaI+gammaI);
x3=(1-sig)/gammaA;

B = [0; 0; R0]; %termini noti
X = [1 0 -betaA_betaP; -betaI_betaA 1 0; x3 x2 x1];
A = linsolve(X,B);
betaA = A(1);
betaI = A(2);
betaP = A(3);
nolockdown= 25; %periodo senza lockdown (supponiamo anche che il cambio di R0 avvenga dopo nolockdown giorni)
novax= 270; % 9 mesi senza vaccinazioni
NV = length(prima_dose_norm); % numero di giorni con vaccini
N= nolockdown + novax + NV;

clear x1 x2 x3 B X A betaA_betaP betaI_betaA NV 


