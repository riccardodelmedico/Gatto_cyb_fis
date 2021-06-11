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
parameters_vaccini;
dati_vaccini;


Lvect = zeros(1,N);
x0= [1-1*E0 1*E0 zeros(1,22)];
options = odeset('RelTol',1e-7,'AbsTol',1e-8);

%cercando il tasso di raddoppio seleziono tradd in t quando ho
%tra gli espostiv E quando raggiungo 2*E0
tic
[t,x_vaccini_tot]= ode45('gatto_vaccini_unico', 0:0.05:10, x0); 
toc


figure(1)
E = x_vaccini_tot(:,2);
plot(t, E);
hold on
plot(t, 2*E0*ones(length(t)));
hold off
legend('Visione raddoppio esposti')
%calcolo tasso di raddoppio t_radd tramite selezioni dei valori di t, con
%il vettore logico di selezione selez

selez = (E< (2*E0 + 1e-8) & E > (2*E0 - 1e-8));
t_radd = t(selez)

%la domanda è: non si
%dovrebbero considerare TUTTI I MALATI? in ogni caso aspettandosi che
%magari il contributi dei restanti compartimenti sia trascurabile o non
%così rilevante?

%vediamo che succede considerando anche i restanti malati, quindi gruppi
%P,I, A, H e Q

figure(2)
malati = E + x_vaccini_tot(:,3) + x_vaccini_tot(:,4) + x_vaccini_tot(:,5) + x_vaccini_tot(:,6) + x_vaccini_tot(:,7);
plot(t, malati);
hold on
plot(t, 2*E0*ones(length(t)));
hold off
legend('Visione raddoppio malati totali')

selez1 = (malati< (2*E0 + 5e-8) & malati > (2*E0 - 5e-8));
t_radd1 = t(selez1)


%% ora cerchiamo di capire che succede se cambio le condizioni iniziali

%l'idea è metetre delle persone in tutti i compartimenti e vedere quando
%raddoppiano

x0 = [1-8*E0, E0*ones(1,8) zeros(1,15) ];

%cercando il tasso di raddoppio seleziono tradd in t quando ho
%tra gli espostiv E quando raggiungo 2*E0
tic
[t,x_vaccini_tot]= ode45('gatto_vaccini_unico', 0:0.02:7, x0,options); 
toc


figure(1)
E = x_vaccini_tot(:,2);
plot(t, E);
hold on
plot(t, 2*E0*ones(length(t)));
hold off
legend('Esposti')

%calcolo tasso di raddoppio t_radd tramite selezioni dei valori di t, con
%il vettore logico di selezione selez

selez = (E< (2*E0 + 1e-8) & E > (2*E0 - 1e-8));
t_radd_E = t(selez)

%e così via vediamo le altre variabili

figure(2)
P = x_vaccini_tot(:,3);
plot(t, P);
hold on
plot(t, 2*E0*ones(length(t)));
hold off
legend('Paucisintomatici')

%calcolo tasso di raddoppio t_radd tramite selezioni dei valori di t, con
%il vettore logico di selezione selez

selez = (P< (2*E0 + 1e-8) & P > (2*E0 - 1e-8));
t_radd_P = t(selez)

%e così via per le altre variabili

figure(3)
I = x_vaccini_tot(:,4);
plot(t, I);
hold on
plot(t, 2*E0*ones(length(t)));
hold off
legend('Infetti')

%calcolo tasso di raddoppio t_radd tramite selezioni dei valori di t, con
%il vettore logico di selezione selez

selez = (I< (2*E0 + 1e-8) & I > (2*E0 - 1e-8));
t_radd_I = t(selez)

%anche se sembra che per raddoppiare ci voglia sempre più tempo

figure(4)
A = x_vaccini_tot(:,5);
plot(t, A);
hold on
plot(t, 2*E0*ones(length(t)));
hold off
legend('Asintomatici')

%calcolo tasso di raddoppio t_radd tramite selezioni dei valori di t, con
%il vettore logico di selezione selez

selez = (A< (2*E0 + 1e-8) & A > (2*E0 - 1e-8));
t_radd_A = t(selez)

%quindi vediamo la somma di tutti i componenti

trasmettitori = P+I+A ;% domanda: ha senso aggiungere anche chi non trasmette?E + x_vaccini_tot(:, 6) + x_vaccini_tot(:, 7);
figure(5)
plot(t, trasmettitori);
hold on
plot(t, 2*E0*3*ones(length(t)));
hold off
legend('Trasmettitori (P,I,A)')

selez = (trasmettitori< (2*3*E0 + 3e-8) & trasmettitori> (2*3*E0 - 3e-8));
t_radd_trasm = t(selez)

figure(6)
malati = E + x_vaccini_tot(:,3) + x_vaccini_tot(:,4) + x_vaccini_tot(:,5) + x_vaccini_tot(:,6) + x_vaccini_tot(:,7);
plot(t, malati);
hold on
plot(t, 2*E0*6*ones(length(t)));
hold off
legend('Malati totali (E,P,I,A,H,Q)')

selez1 = (malati< (2*6*E0 + 5e-8) & malati > (2*6*E0 - 5e-8));
t_radd_tot = t(selez1)


%% rifacciamo l'analisi del tasso di raddoppio nel nuovo modello per vedere se c'è stato un rallentamento 
%causa cambio ipotes i sulle distribuzioni degli infetti (da esponenziale a gamma) 
global x0_casc
x0_casc = [1-1*E0 , 1*E0, zeros(1,29)];

tic
[t,x_vaccini_tot]= ode45('gatto_vaccini_unico_cascate', 0:0.1:10, x0_casc,options); 
toc

%
figure(1)
E = x_vaccini_tot(:,3);
plot(t, E);
hold on
plot(t, 2*E0*ones(length(t)));
hold off
legend('Esposti (distrib. gamma)')

%calcolo tasso di raddoppio t_radd tramite selezioni dei valori di t, con
%il vettore logico di selezione selez

selez = (E< (2*E0 + 1e-8) & E > (2*E0 - 1e-8));
t_radd = t(selez)

%la domanda è: non si
%dovrebbero considerare TUTTI I MALATI? in ogni caso aspettandosi che
%magari il contributi dei restanti compartimenti sia trascurabile o non
%così rilevante?

%vediamo che succede considerando anche i restanti malati, quindi gruppi
 % chi trasmette ossia P,I, A
P = x_vaccini_tot(:,4);
I = x_vaccini_tot(:,5);
A = x_vaccini_tot(:,6);


trasmettitori = P+I+A;
figure(2)
plot(t, trasmettitori);
hold on
plot(t, 2*E0*ones(length(t)));
hold off
legend('Trasmettitori')

selez1 = (trasmettitori< (2*E0 + 1e-8) & trasmettitori > (2*E0 - 1e-8));
t_radd1 = t(selez1);

%qui addirittura sembra essere rilevante solo se consideriamo tutti i
%reparti

figure(3)
malati = E+P+I+A+x_vaccini_tot(:,9)+x_vaccini_tot(:,10);
plot(t, malati);
hold on
plot(t, 2*E0*ones(length(t)));
hold off
legend('Malati totali (E,P,I,A,H,Q)')

selez1 = (malati< (2*E0 + 5e-8) & malati > (2*E0 - 5e-8));
t_radd_tot = t(selez1)
%
%% ora cerchiamo di capire che succede se cambio le condizioni iniziali

%l'idea è metetre delle persone in tutti i compartimenti e vedere quando
%raddoppiano

x0_casc = [1-8*E0, 0, E0, E0, E0 E0  0 0 E0 E0 E0 E0 zeros(1,19) ];

%cercando il tasso di raddoppio seleziono tradd in t quando ho
%tra gli espostiv E quando raggiungo 2*E0
tic
[t,x_vaccini_tot]= ode45('gatto_vaccini_unico_cascate', 0:0.02:7, x0_casc,options); 
toc


figure(1)
E = x_vaccini_tot(:,3);
plot(t, E);
hold on
plot(t, 2*E0*ones(length(t)));
hold off
legend('Esposti')

%calcolo tasso di raddoppio t_radd tramite selezioni dei valori di t, con
%il vettore logico di selezione selez

selez = (E< (2*E0 + 1e-8) & E > (2*E0 - 1e-8));
t_radd_E = t(selez)

%e così via vediamo le altre variabili

figure(2)
P = x_vaccini_tot(:,4);
plot(t, P);
hold on
plot(t, 2*E0*ones(length(t)));
hold off
legend('Paucisintomatici')

%calcolo tasso di raddoppio t_radd tramite selezioni dei valori di t, con
%il vettore logico di selezione selez

selez = (P< (2*E0 + 1e-8) & P > (2*E0 - 1e-8));
t_radd_P = t(selez)

%e così via per le altre variabili

figure(3)
I = x_vaccini_tot(:,5);
plot(t, I);
hold on
plot(t, 2*E0*ones(length(t)));
hold off
legend('Infetti')

%calcolo tasso di raddoppio t_radd tramite selezioni dei valori di t, con
%il vettore logico di selezione selez

selez = (I< (2*E0 + 1e-8) & I > (2*E0 - 1e-8));
t_radd_I = t(selez)

%anche se sembra che per raddoppiare ci voglia sempre più tempo

figure(4)
A = x_vaccini_tot(:,6);
plot(t, A);
hold on
plot(t, 2*E0*ones(length(t)));
hold off
legend('Asintomatici')

%calcolo tasso di raddoppio t_radd tramite selezioni dei valori di t, con
%il vettore logico di selezione selez

selez = (A< (2*E0 + 1e-8) & A > (2*E0 - 1e-8));
t_radd_A = t(selez)

%quindi vediamo la somma di tutti i componenti

trasmettitori = P+I+A ;% domanda: ha senso aggiungere anche chi non trasmette?E + x_vaccini_tot(:, 6) + x_vaccini_tot(:, 7);
figure(5)
plot(t, trasmettitori);
hold on
plot(t, 2*E0*3*ones(length(t)));
hold off
legend('Trasmettitori (P,I,A)')

selez = (trasmettitori< (2*3*E0 + 3e-8) & trasmettitori> (2*3*E0 - 3e-8));
t_radd_trasm = t(selez)

figure(6)
malati = E + P + I + A + x_vaccini_tot(:, 9) + x_vaccini_tot(:, 10);
plot(t, malati);
hold on
plot(t, 2*E0*6*ones(length(t)));
hold off
legend('Malati totali (E,P,I,A,H,Q)')

selez1 = (malati< (2*6*E0 + 5e-8) & malati > (2*6*E0 - 5e-8));
t_radd_tot = t(selez1)



%% rifacciamo l'analisi del tasso di raddoppio nel nuovo modello per vedere se c'è stato un rallentamento 
%causa cambio ipotes i sulle distribuzioni degli infetti (da esponenziale a gamma) 
global x0_casc_inf
x0_casc_inf = [1-1*E0 , 1*E0, zeros(1,28)];

tic
[t,x_vaccini_tot]= ode45('gatto_vaccini_unico_cascatesoloInfetti', 0:0.1:10, x0_casc_inf,options); 
toc

%
figure(1)
E = x_vaccini_tot(:,2);
plot(t, E);
hold on
plot(t, 2*E0*ones(length(t)));
hold off
legend('Esposti')

%calcolo tasso di raddoppio t_radd tramite selezioni dei valori di t, con
%il vettore logico di selezione selez

selez = (E< (2*E0 + 1e-8) & E > (2*E0 - 1e-8));
t_radd = t(selez)

%la domanda è: non si
%dovrebbero considerare TUTTI I MALATI? in ogni caso aspettandosi che
%magari il contributi dei restanti compartimenti sia trascurabile o non
%così rilevante?

%vediamo che succede considerando anche i restanti malati, quindi gruppi
%P,I, A, H e Q

figure(2)
malati = E + x_vaccini_tot(:,3) + x_vaccini_tot(:,6) + x_vaccini_tot(:,7) +...
         x_vaccini_tot(:,8)+ x_vaccini_tot(:,9);
plot(t, malati);
hold on
plot(t, 2*E0*ones(length(t)));
hold off
legend('Malati')

selez1 = (malati< (2*E0 + 1e-8) & malati > (2*E0 - 1e-8));
t_radd1 = t(selez1)

figure(3)
trasmettitori = x_vaccini_tot(:,3) +  x_vaccini_tot(:,6) +  x_vaccini_tot(:,7);
plot(t, trasmettitori);
hold on
plot(t, 2*E0*ones(length(t)));
hold off
legend('Trasmettitori')

selez1 = (trasmettitori< (2*E0 + 5e-8) & trasmettitori > (2*E0 - 5e-8));
t_radd_tot = t(selez1)

%



%% ora cerchiamo di capire che succede se cambio le condizioni iniziali

%l'idea è metetre delle persone in tutti i compartimenti e vedere quando
%raddoppiano

x0_casc_inf = [1-7*E0, E0, E0, 0 0 E0 E0 E0 E0 E0 zeros(1,20) ];

%cercando il tasso di raddoppio seleziono tradd in t quando ho
%tra gli espostiv E quando raggiungo 2*E0
tic
[t,x_vaccini_tot]= ode45('gatto_vaccini_unico_cascatesoloInfetti', 0:0.02:7, x0_casc_inf,options); 
toc


figure(1)
E = x_vaccini_tot(:,2);
plot(t, E);
hold on
plot(t, 2*E0*ones(length(t)));
hold off
legend('Esposti')

%calcolo tasso di raddoppio t_radd tramite selezioni dei valori di t, con
%il vettore logico di selezione selez

selez = (E< (2*E0 + 1e-8) & E > (2*E0 - 1e-8));
t_radd_E = t(selez)

%e così via vediamo le altre variabili

figure(2)
P = x_vaccini_tot(:,3);
plot(t, P);
hold on
plot(t, 2*E0*ones(length(t)));
hold off
legend('Paucisintomatici')

%calcolo tasso di raddoppio t_radd tramite selezioni dei valori di t, con
%il vettore logico di selezione selez

selez = (P< (2*E0 + 1e-8) & P > (2*E0 - 1e-8));
t_radd_P = t(selez)

%e così via per le altre variabili

figure(3)
I = x_vaccini_tot(:,6);
plot(t, I);
hold on
plot(t, 2*E0*ones(length(t)));
hold off
legend('Infetti')

%calcolo tasso di raddoppio t_radd tramite selezioni dei valori di t, con
%il vettore logico di selezione selez

selez = (I< (2*E0 + 1e-8) & I > (2*E0 - 1e-8));
t_radd_I = t(selez)

%anche se sembra che per raddoppiare ci voglia sempre più tempo

figure(4)
A = x_vaccini_tot(:,7);
plot(t, A);
hold on
plot(t, 2*E0*ones(length(t)));
hold off
legend('Asintomatici')

%calcolo tasso di raddoppio t_radd tramite selezioni dei valori di t, con
%il vettore logico di selezione selez

selez = (A< (2*E0 + 1e-8) & A > (2*E0 - 1e-8));
t_radd_A = t(selez)

%quindi vediamo la somma di tutti i componenti

trasmettitori = P+I+A ;% domanda: ha senso aggiungere anche chi non trasmette?E + x_vaccini_tot(:, 6) + x_vaccini_tot(:, 7);
figure(5)
plot(t, trasmettitori);
hold on
plot(t, 2*E0*3*ones(length(t)));
hold off
legend('Trasmettitori (P,I,A)')

selez = (trasmettitori< (2*3*E0 + 3e-8) & trasmettitori> (2*3*E0 - 3e-8));
t_radd_trasm = t(selez)

figure(6)
malati = E + P + I + A + x_vaccini_tot(:, 9) + x_vaccini_tot(:, 8);
plot(t, malati);
hold on
plot(t, 2*E0*6*ones(length(t)));
hold off
legend('Malati totali (E,P,I,A,H,Q)')

selez1 = (malati< (2*6*E0 + 5e-8) & malati > (2*6*E0 - 5e-8));
t_radd_tot = t(selez1)
