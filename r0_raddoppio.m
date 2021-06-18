%% Script per il cambio di R0, da settare nel main
global lambda deltaE deltaP sig eta gammaI alfaI gammaA zeta gammaH ...
    alfaH gammaQ gammaA betaP betaI betaA eff1 eff2 ef1 teta R0
%i seguenti per il modello base
%R0 = 2.13 per raddoppio pari a circa 3.5
%R0 = 2.6 per raddoppio a circa 2.5
%R0 = 2.3 per una via di mezzo, con tempo raddoppio pari a 3

%i seguenti per il modello realistico del gatto (cascate E,H)
%R0 = 2.3 per tempo raddoppio pari a 3.5
%R0 = 2.6  per tempo raddoppio pari a 3
%R0 = 2.95 per tempo raddoppio pari a 2.5

%per il modello realistico di Manfredi (cascate E,I,A) 
%R0 = 2.58 per raddoppio pari a circa 3.5
%R0 = 3.3 per raddoppio a circa 2.5
%R0 = 2.9 per una via di mezzo, con tempo raddoppio pari a 3

R0=1;

betaA_betaP= 0.033;
betaI_betaA= 1.03;
x1=1/deltaP;
x2=sig/(eta+alfaI+gammaI);
x3=(1-sig)/gammaA;

B = [0; 0; R0]; %termini noti
X = [1 0 -betaA_betaP; -betaI_betaA 1 0; x3 x2 x1];
A = linsolve(X,B);
betaA = A(1);
betaI = A(2);
betaP = A(3);

clear x1 x2 x3 B X A betaA_betaP betaI_betaA