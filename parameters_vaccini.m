global lambda deltaE deltaP sigma eta gammaI alfaI gammaA zeta gammaH ...
       alfaH gammaQ gammaA betaP betaI betaA eff1 eff2 ef1 theta

R0=3.6;
deltaE= 1 / 3.32;
deltaP= 1 / 0.75;
eta= 1 / 4.05;
gammaI= 1 / 14.32;
gammaQ=gammaI;
gammaH= gammaI;
gammaA= 2*gammaI;
alfaI= 1 / 24.23;
sigma= 0.25; % (1-sigma) frazione degli asintomatici
betaP = R0/(1/deltaP + 1.03*sigma/( eta + alfaI + gammaI)+ 0.033*(1-sigma) / gammaA);
betaA= 0.033*betaP; %betaA asymptomatic transmission rate
betaI= 1.03*betaA;
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
theta = 0.5; %efficacia del lockdown
% options = optimoptions('fmincon','Display','none','Algorithm','active-set',...
%     'OptimalityTolerance', 1e-1, 'MaxFunctionEvaluations', 5,'FunValCheck','on');