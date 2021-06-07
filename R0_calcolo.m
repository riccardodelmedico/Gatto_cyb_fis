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


%% sperimentale

%%

%cercando il tasso di raddoppio seleziono tradd in t quando ho
%tra gli espostiv E quando raggiungo 2*E0
tic
[t,x_vaccini_tot]= ode45('gatto_vaccini_unico', 0:0.1:5, x0,options); 
toc


figure(3)
E = x_vaccini_tot(:,2);
plot(t, E);
hold on
plot(t, 2*E0*ones(length(t)));
hold off

%calcolo tasso di raddoppio t_radd tramite selezioni dei valori di t, con
%il vettore logico di selezione selez

selez = (E< (2*E0 + 1e-8) & E > (2*E0 - 1e-8));
t_radd = t(selez)

%la domanda �: non si
%dovrebbero considerare TUTTI I MALATI? in ogni caso aspettandosi che
%magari il contributi dei restanti compartimenti sia trascurabile o non
%cos� rilevante?

%vediamo che succede considerando anche i restanti malati, quindi gruppi
%P,I, A, H e Q

figure(4)
malati = E + x_vaccini_tot(:,3) + x_vaccini_tot(:,4) + x_vaccini_tot(:,5) + x_vaccini_tot(:,6) + x_vaccini_tot(:,7);
plot(t, malati);
hold on
plot(t, 2*E0*ones(length(t)));
hold off

selez1 = (malati< (2*E0 + 1e-8) & malati > (2*E0 - 1e-8));
t_radd1 = t(selez1)


%% ora cerchiamo di capire che succede se cambio le condizioni iniziali

%l'idea � metetre delle persone in tutti i compartimenti e vedere quando
%raddoppiano

x0 = [1-E0, E0*ones(1,8) zeros(1,15) ];

%cercando il tasso di raddoppio seleziono tradd in t quando ho
%tra gli espostiv E quando raggiungo 2*E0
tic
[t,x_vaccini_tot]= ode45('gatto_vaccini_unico', 0:0.1:5, x0,options); 
toc


figure(1)
E = x_vaccini_tot(:,2);
plot(t, E);
hold on
plot(t, 2*E0*ones(length(t)));
hold off

%calcolo tasso di raddoppio t_radd tramite selezioni dei valori di t, con
%il vettore logico di selezione selez

selez = (E< (2*E0 + 1e-8) & E > (2*E0 - 1e-8));
t_radd_E = t(selez)

%e cos� via vediamo le altre variabili

figure(2)
P = x_vaccini_tot(:,3);
plot(t, P);
hold on
plot(t, 2*E0*ones(length(t)));
hold off

%calcolo tasso di raddoppio t_radd tramite selezioni dei valori di t, con
%il vettore logico di selezione selez

selez = (P< (2*E0 + 1e-8) & P > (2*E0 - 1e-8));
t_radd_P = t(selez)

%e cos� via per le altre variabili

figure(3)
I = x_vaccini_tot(:,4);
plot(t, I);
hold on
plot(t, 2*E0*ones(length(t)));
hold off

%calcolo tasso di raddoppio t_radd tramite selezioni dei valori di t, con
%il vettore logico di selezione selez

selez = (I< (2*E0 + 1e-8) & I > (2*E0 - 1e-8));
t_radd_I = t(selez)

%anche se sembra che per raddoppiare ci voglia sempre pi� tempo

figure(4)
A = x_vaccini_tot(:,5);
plot(t, A);
hold on
plot(t, 2*E0*ones(length(t)));
hold off

%calcolo tasso di raddoppio t_radd tramite selezioni dei valori di t, con
%il vettore logico di selezione selez

selez = (A< (2*E0 + 1e-8) & A > (2*E0 - 1e-8));
t_radd_A = t(selez)

%quindi vediamo la somma di tutti i componenti

malati = E+P+I+A ;% domanda: ha senso aggiungere anche chi non trasmette?+ x_vaccini_tot(:, 6) + x_vaccini_tot(:, 7);
figure(5)
plot(t, malati);
hold on
plot(t, 2*E0*4*ones(length(t)));
hold off

selez = (malati< (2*4*E0 + 3e-8) & malati> (2*4*E0 - 3e-8));
t_radd_tot = t(selez)


%%%% rifacciamo l'analisi del tasso di raddoppio nel nuovo modello per vedere se c'� stato un rallentamento 
%causa cambio ipotes i sulle distribuzioni degli infetti (da esponenziale a gamma) 

tic
[t,x_vaccini_tot]= ode45('gatto_vaccini_unico_cascate', 0:0.1:10, x0_casc,options); 
toc

%
figure(3)
E = x_vaccini_tot(:,3);
plot(t, E);
hold on
plot(t, 2*E0*ones(length(t)));
hold off

%calcolo tasso di raddoppio t_radd tramite selezioni dei valori di t, con
%il vettore logico di selezione selez

selez = (E< (2*E0 + 1e-8) & E > (2*E0 - 1e-8));
t_radd = t(selez)

%la domanda �: non si
%dovrebbero considerare TUTTI I MALATI? in ogni caso aspettandosi che
%magari il contributi dei restanti compartimenti sia trascurabile o non
%cos� rilevante?

%vediamo che succede considerando anche i restanti malati, quindi gruppi
%P,I, A, H e Q

figure(4)
malati = E + x_vaccini_tot(:,2) + x_vaccini_tot(:,4) + x_vaccini_tot(:,5) + x_vaccini_tot(:,6) + x_vaccini_tot(:,7) +...
         x_vaccini_tot(:,8)+ x_vaccini_tot(:,9) +x_vaccini_tot(:,10) + x_vaccini_tot(:,11) + x_vaccini_tot(:,12) ;
plot(t, malati);
hold on
plot(t, 2*E0*ones(length(t)));
hold off

selez1 = (malati< (2*E0 + 1e-8) & malati > (2*E0 - 1e-8));
t_radd1 = t(selez1);

%qui addirittura sembra essere rilevante solo se consideriamo tutti i
%reparti

%%
%% rifacciamo l'analisi del tasso di raddoppio nel nuovo modello per vedere se c'� stato un rallentamento 
%causa cambio ipotes i sulle distribuzioni degli infetti (da esponenziale a gamma) 

tic
[t,x_vaccini_tot]= ode45('gatto_vaccini_unico_cascatesoloInfetti', 0:0.1:10, x0_casc_inf,options); 
toc

%
figure(3)
E = x_vaccini_tot(:,3);
plot(t, E);
hold on
plot(t, 2*E0*ones(length(t)));
hold off

%calcolo tasso di raddoppio t_radd tramite selezioni dei valori di t, con
%il vettore logico di selezione selez

selez = (E< (2*E0 + 1e-8) & E > (2*E0 - 1e-8));
t_radd = t(selez)

%la domanda �: non si
%dovrebbero considerare TUTTI I MALATI? in ogni caso aspettandosi che
%magari il contributi dei restanti compartimenti sia trascurabile o non
%cos� rilevante?

%vediamo che succede considerando anche i restanti malati, quindi gruppi
%P,I, A, H e Q

figure(4)
malati = E + x_vaccini_tot(:,2) + x_vaccini_tot(:,4) + x_vaccini_tot(:,5) + x_vaccini_tot(:,6) + x_vaccini_tot(:,7) +...
         x_vaccini_tot(:,8)+ x_vaccini_tot(:,9) +x_vaccini_tot(:,10) + x_vaccini_tot(:,11) ;
plot(t, malati);
hold on
plot(t, 2*E0*ones(length(t)));
hold off

selez1 = (malati< (2*E0 + 1e-8) & malati > (2*E0 - 1e-8));
t_radd1 = t(selez1)
