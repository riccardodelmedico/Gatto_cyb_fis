%%R0 NGM
parameters_vaccini;

J0=zeros(4,4);
T=zeros(4,4);
sigma_maius=zeros(4,4);

J0=[-deltaE, betaP, betaI, betaA; deltaE, -deltaP, 0, 0; 0, sigma*deltaP, ...
    -(eta+alfaI+gammaI), 0; 0, (1-sigma)*deltaP, 0, -gammaA];

T=[0, betaP, betaI, betaA; 0,0,0,0;0,0,0,0;0,0,0,0];

sigma_maius= [-deltaE, 0,0,0; deltaE, -deltaP, 0, 0; 0, sigma*deltaP, ...
    -(eta+alfaI+gammaI), 0; 0 (1-sigma)*deltaP, 0, -gammaA];

NGM=-T*inv(sigma_maius);

R0=max(abs(eig(NGM)))

%% R0 
%scrivo un equilibrio endemico

%syms lambda deltaE deltaP sigma eta gammaI alfaI gammaA zeta alfaH gammaH gammaQ
format longg
A=zeros(7,7);
A(1,1)= -lambda;
A(2,1)= lambda;
A(2,2)= -deltaE;
A(3,3)= deltaE;
A(3,4)= - deltaP;
A(4,3)= sigma*deltaP;
A(4,4)= - (eta + gammaI + alfaI);
A(5,3)= (1-sigma)*deltaP;
A(5,5)= gammaA;
A(6,4)= (1-zeta)*eta;
A(6,6)= -(gammaH + alfaH);
A(7,4)= zeta*eta;
A(7,7)= -gammaQ;
% A(8,4)= gammaI;
% A(8,5)= gammaA;
% A(8,6)= gammaH;
% A(9,4)= alfaI;
% A(9,6)= alfaH;

B=[1;0;0;0;0;0;0];
sol= linsolve(A,B);

A*sol

