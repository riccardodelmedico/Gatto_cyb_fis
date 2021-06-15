global lambda deltaE deltaP sig eta gammaI alfaI gammaA zeta gammaH ...
    alfaH gammaQ gammaA betaP betaI betaA eff1 eff2 ef1 teta

%R0 = 2.13 per raddoppio pari a circa 3.5
%R0 = 2.6 per raddoppio a circa 2.5
%R0 = 2.3 per una via di mezzo, con tempo raddoppio pari a 3

R0=2.3;
deltaE= 1 / 3.32;
deltaP= 1 / 0.75;
eta= 1 / 4.05;
gammaI= 1 / 14.32;
gammaQ=gammaI;
gammaH= gammaI;
gammaA= 2*gammaI;
alfaI= 1 / 24.23;
sig= 0.25; % (1-sigma) frazione degli asintomatici
%betaP = R0/(1/deltaP + 1.03*sig/( eta + alfaI + gammaI)+ 0.033*(1-sig) / gammaA);
betaA_betaP= 0.033; %betaA asymptomatic transmission rate
betaI_betaA= 1.03;
% betaP1/ betaP== 0.82; % betaP1 after restriction on february 22, 2020
% betaP2/ betaP1== 0.66; % betaP2 after restriction on march 8, 2020
% deltat0= 34.94;
% omega= 7.84;
zeta= 0.4; % fraction of sympomatic infected being quarantined
alfaH = alfaI;
% rs= 0.5;
eff1 = 0.9;
eff2 = 1;
ef1 = 0.7;
teta = 1; %efficacia del lockdown
% options = optimoptions('fmincon','Display','none','Algorithm','active-set',...
%     'OptimalityTolerance', 1e-1, 'MaxFunctionEvaluations', 5,'FunValCheck','on');

x1=1/deltaP;
x2=sig/(eta+alfaI+gammaI);
x3=(1-sig)/gammaA;

B = [0 0 R0].'; %termini noti
X = [1 0 -betaA_betaP; -betaI_betaA 1 0; x3 x2 x1];
A = linsolve(X,B);
betaA = A(1);
betaI = A(2);
betaP = A(3);
nolockdown= 25;
novax= 270;
Mvax = size(prima_dose_norm);
NV = Mvax(1);
N= nolockdown + novax + NV;