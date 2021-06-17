%%R0 NGM
global lambda deltaE deltaP sig eta gammaI alfaI gammaA zeta gammaH ...
       alfaH gammaQ gammaA x0 N eff1 eff2 ef1 prima_d seconda_d Lvect NV
parameters_vaccini;
dati_vaccini;

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


%% sperimentale

global lambda deltaE deltaP sig eta gammaI alfaI gammaA zeta gammaH ...
       alfaH gammaQ gammaA x0 N eff1 eff2 ef1 prima_d seconda_d Lvect NV

format longg
parameters_vaccini;
dati_vaccini;

Lvect = zeros(1,N);
x0= [1-1*E0 1*E0 zeros(1,22)];
options = odeset('RelTol',1e-9,'AbsTol',1e-9);

step = 0.1;
%cercando il tasso di raddoppio seleziono tradd in t quando ho
%tra gli espostiv E quando raggiungo 2*E0
tic
[t,x_vaccini_tot]= ode45('gatto_vaccini_unico', 0:step:10, x0);
x_R0 = x_vaccini_tot(end, :);
[t1,x_vaccini_tot1]= ode45('gatto_vaccini_unico', 0:step:20, x_R0);
toc
% tic
% [yout]= ode4('gatto_vaccini_unico', 0, 0.1, 30, x0)
% toc
%VEDERE CHE CAZZO VUOLE E PERCHè DA I DATI NEGATIVI, QUESTA SEZIONE è
%L'UNICA PROBLEMATICA PER QUESTO FATTO ED è L'UNICA CHE FA CASINO
%ALL'INIZIO

t = 0:step:30+step;

E = [x_vaccini_tot(:,2); x_vaccini_tot1(:,2)];
P = [x_vaccini_tot(:,3); x_vaccini_tot1(:,3)];
I = [x_vaccini_tot(:,4);x_vaccini_tot1(:,4)];
A = [x_vaccini_tot(:,5);x_vaccini_tot1(:,5)];

figure(1)
plot(t, log(E))
hold on
plot(t,log(P))
plot(t,log(I))
plot(t,log(A))
hold off
legend('E', 'P','I', 'A')
xlabel('Days')
ylabel('Logarithmic scale')

figure(2)
plot(t,[E P I A])
legend('E', 'P','I', 'A')
xlabel('Days')


%calcolo il tasso esponenziale di crescita delle curve come coefficiente angolare delle curve logaritmiche a regime, 
%che poi uso per il tempo di raddoppio--->problema: come collego Rt al
%coefficiente angolare? NON SI PUò FARE ANALITICAMENTE

tasso_expo = (log(P(end)) - log(P(end-1)))/step
t_radd = log(2)/tasso_expo


%% ritaratura R0, mettendo il tempo di raddoppio a CIRCA 3, A ''OCCHIO''

parameters_vaccini_R0_raddoppio;
dati_vaccini;


Lvect = zeros(1,N);
x0= [1-1*E0 1*E0 zeros(1,22)];
options = odeset('RelTol',1e-9,'AbsTol',1e-9);

step = 0.1;
%cercando il tasso di raddoppio seleziono tradd in t quando ho
%tra gli espostiv E quando raggiungo 2*E0
tic
[t,x_vaccini_tot]= ode45('gatto_vaccini_unico', 0:step:10, x0);
x_R0 = x_vaccini_tot(end, :);
[t1,x_vaccini_tot1]= ode45('gatto_vaccini_unico', 0:step:20, x_R0);
toc

%VEDERE CHE CAZZO VUOLE E PERCHè DA I DATI NEGATIVI, QUESTA SEZIONE è
%L'UNICA PROBLEMATICA PER QUESTO FATTO ED è L'UNICA CHE FA CASINO
%ALL'INIZIO

t = 0:step:30+step;

E = [x_vaccini_tot(:,2); x_vaccini_tot1(:,2)];
P = [x_vaccini_tot(:,3); x_vaccini_tot1(:,3)];
I = [x_vaccini_tot(:,4);x_vaccini_tot1(:,4)];
A = [x_vaccini_tot(:,5);x_vaccini_tot1(:,5)];

figure(1)
plot(t, log(E))
hold on
plot(t,log(P))
plot(t,log(I))
plot(t,log(A))
hold off
legend('E', 'P','I', 'A')

figure(2)
plot(t,[E P I A])

%calcolo il tasso esponenziale di crescita delle curve come coefficiente angolare delle curve logaritmiche a regime, 
%che poi uso per il tempo di raddoppio--->problema: come collego Rt al
%coefficiente angolare? NON SI PUò FARE ANALITICAMENTE

tasso_expo = (log(P(end)) - log(P(end-1)))/step
t_radd = log(2)/tasso_expo

%OSS: NON SERVE CAMBIARE LE CONDIZIONI INIZIALI SE CI INTERESSA
%L'EVOLUZIONE ESPONENZIALE DOPO UN PO' DI TEMPO (UNA SPECIE DI REGIME INIZIALE)
% %% ora cerchiamo di capire che succede se cambio le condizioni iniziali
% 
% %l'idea è metetre delle persone in tutti i compartimenti e vedere quando
% %raddoppiano
% 
% x0 = [1-E0, E0*ones(1,8) zeros(1,15) ];
% 
% step = 0.5;
% %cercando il tasso di raddoppio seleziono tradd in t quando ho
% %tra gli espostiv E quando raggiungo 2*E0
% tic
% [t,x_vaccini_tot]= ode45('gatto_vaccini_unico', 0:step:30, x0,options); 
% toc
% 
% 
% 
% % facciamo i plot logaritmici tanto dobbiamo vedere il funzionamento a
% % regime dell'epidemia, in maniera esponenziale
% E = x_vaccini_tot(:,2);
% P = x_vaccini_tot(:,3);
% I = x_vaccini_tot(:,4);
% A = x_vaccini_tot(:,5);
% 
% 
% figure(3)
% plot(t, log(E))
% hold on
% plot(t,log(P))
% plot(t,log(I))
% plot(t,log(A))
% hold off
% legend('E', 'P','I', 'A')
% 
% figure(4)
% plot(t,[E P I A])
% 
% tasso_expo = (log(E(end)) - log(E(end-1)))/step
% t_radd = log(2)/tasso_expo



%% rifacciamo l'analisi del tasso di raddoppio nel nuovo modello per vedere se c'è stato un rallentamento 
%causa cambio ipotes i sulle distribuzioni degli infetti (da esponenziale a gamma) 

%rifacciamo lo studio originale

parameters_vaccini
dati_vaccini


global x0_casc
x0_casc = [1-1*E0 , 1*E0, zeros(1,29)];
step = 0.1;

tic
[t,x_vaccini_tot]= ode45('gatto_vaccini_unico_cascate', 0:step:25, x0_casc,options); 
toc
%ricorda che qui su E abbiamo la distribuzione gamma a due compartimenti

E = x_vaccini_tot(:,3);
P = x_vaccini_tot(:,4);
I = x_vaccini_tot(:,5);
A = x_vaccini_tot(:,6);


figure(5)
plot(t, log(E))
hold on
plot(t,log(P))
plot(t,log(I))
plot(t,log(A))
hold off
legend('E', 'P','I', 'A')

figure(6)
plot(t,[E P I A])

tasso_expo = (log(P(end)) - log(P(end-1)))/step
t_radd = log(2)/tasso_expo

%% cerchiamo di tarare allora 
parameters_vaccini_R0_raddoppio
dati_vaccini


global x0_casc
x0_casc = [1-1*E0 , 1*E0, zeros(1,29)];
step = 0.1;

tic
[t,x_vaccini_tot]= ode45('gatto_vaccini_unico_cascate', 0:step:25, x0_casc,options); 
toc
%ricorda che qui su E abbiamo la distribuzione gamma a due compartimenti

E = x_vaccini_tot(:,3);
P = x_vaccini_tot(:,4);
I = x_vaccini_tot(:,5);
A = x_vaccini_tot(:,6);


figure(5)
plot(t, log(E))
hold on
plot(t,log(P))
plot(t,log(I))
plot(t,log(A))
hold off
legend('E', 'P','I', 'A')

figure(6)
plot(t,[E P I A])

tasso_expo = (log(P(end)) - log(P(end-1)))/step
t_radd = log(2)/tasso_expo



%% rifacciamo l'analisi del tasso di raddoppio nel nuovo modello per vedere se c'è stato un rallentamento 
%causa cambio ipotes i sulle distribuzioni degli infetti (da esponenziale a gamma) 
parameters_vaccini

global x0_casc_inf
x0_casc_inf = [1-1*E0 , 1*E0, zeros(1,40)];

step = 0.3;
tic
[t,x_vaccini_tot]= ode45('gatto_vaccini_unico_cascatesoloInfetti', 0:step:25, x0_casc_inf,options); 
toc

E = x_vaccini_tot(:,2);
P = x_vaccini_tot(:,3);
I = x_vaccini_tot(:,6);
A = x_vaccini_tot(:,7);


figure(7)
plot(t, log(E))
hold on
plot(t,log(P))
plot(t,log(I))
plot(t,log(A))
hold off
legend('E', 'P','I', 'A')


figure(8)
plot(t,[E P I A])

tasso_expo = (log(P(end)) - log(P(end-1)))/step
t_radd = log(2)/tasso_expo

%% ritariamo anche questo

parameters_vaccini_R0_raddoppio

global x0_casc_inf
x0_casc_inf = [1-1*E0 , 1*E0, zeros(1,40)];

step = 0.3;
tic
[t,x_vaccini_tot]= ode45('gatto_vaccini_unico_cascatesoloInfetti', 0:step:25, x0_casc_inf,options); 
toc

E = x_vaccini_tot(:,2);
P = x_vaccini_tot(:,3);
I = x_vaccini_tot(:,6);
A = x_vaccini_tot(:,7);


figure(7)
plot(t, log(E))
hold on
plot(t,log(P))
plot(t,log(I))
plot(t,log(A))
hold off
legend('E', 'P','I', 'A')


figure(8)
plot(t,[E P I A])

tasso_expo = (log(P(end)) - log(P(end-1)))/step
t_radd = log(2)/tasso_expo

