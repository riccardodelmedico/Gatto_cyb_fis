clc
clear
clear all
%%
global lambda deltaE deltaP sig eta gammaI alfaI gammaA zeta gammaH ...
       alfaH gammaQ gammaA x0 N eff1 eff2 ef1 prima_d seconda_d Lvect NV

   
format longg
parameters_vaccini;
dati_vaccini;


% prima_d=prima_dose_norm;
% seconda_d=seconda_dose_norm;

% Caricamento dati iniziali opzione 1
%si fa evolvere l'epidemia con controllo dopo 18 giorni(lockdown)
%e senza usare vaccini per un numero specificato di giorni, poi si comincia
%l'intervento con i vaccini con le somministrazioni giornaliere effettuate
%in Italia a partire dal 27/12/2020
nolockdown= 25;
novax= 270;
N= nolockdown + novax + NV;
time= 0:1:N-1;
E0= 10/pop;
 
prima_d= [zeros(novax+nolockdown,1); prima_dose_norm];
seconda_d=[zeros(novax+nolockdown,1); seconda_dose_norm];
%
figure(1)
plot(time, prima_d)
hold on
plot(time,seconda_d)
hold off

%
% %% Caricamento dati iniziali opzione 2
% %si utilizzano i valori reali dell'epidemia nel giorno in cui sono
% %cominciate le vaccinazioni (27/12/2020) (dati da
% %https://github.com/pcm-dpc/COVID-19/blob/master/dati-andamento-nazionale/dpc-covid19-ita-andamento-nazionale-20201227.csv)
% %ipotesi: gli ospedalizzati comprendo TI e ricoverati
% %ipotesi: gli attualmente positivi sono divisi al 30% asintomatici,
% % 20% infettivi e 50% presintomatici
% N=136;
% time= 0:1:N-1;
% totale_ospedalizzati= 26151; ( 0.000441310734902611 per la popolazione normalizzata)
% totale_quarantena= 555609; (0.00937616978733146 idem)
% totale_positivi= 581760;  (0.00981748052223407 idem)
% deceduti= 71925;  (0.00121376905693359 idem)
% dimessi_guariti= 1394011; %sarebbero i nostri recuperati (0.0235246078112624 idem)


%%
%VEDIAMO CHE SUCCEDE ACCOPPIANDO LE DINAMICHE
parameters_vaccini %i parametri sono quelli originali
options = odeset('RelTol',1e-7,'AbsTol',1e-8);
Lvect= [zeros(nolockdown,1)];% 0.5*ones(novax,1); 0.5*ones(136,1)]; %si suppone che il lockdown
%duri a partire dal 25esimo giorno e termini il giorno in cui cominciano le
%vaccinazioni
x0= [1-1*E0 1*E0 zeros(1,22)];

%risolviamo allora solo sull'intervallo di nolockdown, per cercare il tasso
%di raddoppio dell'epidemia libera

tic
[t,x_vaccini_tot]= ode45('gatto_vaccini_unico', 0:1:nolockdown-1, x0,options); 
toc

figure(1)
plot(t, x_vaccini_tot(:,1)) %plot totale

% figure(2)
% plot(t, x_vaccini_tot(:,[1 10 19])) %visione suscettibili
% legend('S(t)','S1(t)','S2(t)')
%
figure(2)
plot(t, x_vaccini_tot(:,2:9)) %altre variabili 1 gatto
legend('E(t)', 'P(t)', 'I(t)', 'A(t)', 'H(t)', 'Q(t)', 'R(t)', 'D(t)')
%
% figure(4) %altre variabili 2 gatto
% plot(t, x_vaccini_tot(:,11:18))
% legend('E1(t)', 'P1(t)', 'I1(t)', 'A1(t)', 'H1(t)', 'Q1(t)', 'R1(t)', 'D1(t)')
% 
% figure(5) %altre variabili 3 gatto
% plot(t, x_vaccini_tot(:,20:24))
% legend('S2(t)','E2(t)', 'P2(t)', 'I2(t)', 'A2(t)','R2(t)')
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


%% ora cambiamo R0, cos� scalando tutto e vediamo come evolve il sistema
% carichiamo i nuovi parametri dei vaccini 
%(si potrebbe fare anche con il lockdown, ma cos� scaliamo direttamente R0 e i beta)
tic
[t,x_vaccini_tot]= ode45('gatto_vaccini_unico', 0:1:nolockdown-1, x0,options); 
toc

parameters_vaccini_R0_1;
Lvect = zeros(N,1);

%aggiorniamo le condizioni iniziali del nuovo sistema
x0 = x_vaccini_tot(end,:);% variazione suscettibili

%risolviamo il Gatto senza vaccini per 6 mesi
tic
[t1,x_vaccini_tot1]= ode45('gatto_vaccini_unico', 0:1:novax, x0,options); 
toc

figure (1)
t_tot = 0:1:(novax+nolockdown);
plot(t_tot, [x_vaccini_tot(:,1); x_vaccini_tot1(:,1)]) %visione suscettibili

figure(2)
plot(t_tot, [x_vaccini_tot(:,2:end); x_vaccini_tot1(:,2:end)]) %visione restanti



 %% equivalente scalatura di R0, assunto che il lockdown abbia contributo quadratico
 %semplicemente vogliamo una scalatura di 3.6 del termine lambda, tramite
 %la moltiplicazione di (1-L)^2 = 1/3.6=0.2777---->(1-L)= 0.527--->L = 0.473
 
%rimettiamo allora R0 a 3.6
parameters_vaccini 
Lvect = [zeros(nolockdown,1); 0.473*ones(N-nolockdown,1) ];

%rimettiamo le condizioni iniziali come prima
x0= [1-1*E0 1*E0 zeros(1,22)];

%risolviamo il Gatto senza vaccini per tutti i mesi prima dei vaccini
tic
[t1,x_vaccini_tot2]= ode45('gatto_vaccini_unico', 0:1:(novax+nolockdown-1), x0,options); 
toc

figure (1)
t_tot = 0:1:(novax+nolockdown-1);
plot(t_tot, x_vaccini_tot2(:,1))  % visione suscettibili

figure(2)
plot(t_tot, x_vaccini_tot2(:,2:end))  % visione restanti variabili gatto


%% ora introduciamo le cascate di Ode ma secondo il Gatto (su E e su H come a pgina 11)
global x0_casc

parameters_vaccini;

Lvect= [zeros(nolockdown,1); 0.4*ones(novax,1); 0.4*ones(136,1)]; %si suppone che il lockdown
%duri a partire dal 18esimo giorno e termini il giorno in cui cominciano le
%vaccinazioni
x0_casc = [1-1*E0 , 1*E0, zeros(1,29)];
options = odeset('RelTol',1e-8,'AbsTol',1e-10);
tic
[t,x_vaccini_tot]= ode45('gatto_vaccini_unico_cascate', 0:1:nolockdown+novax-1, x0_casc,options); 
toc

figure(1)
plot(t, x_vaccini_tot(:,:)) %plot totale

figure(2)
plot(t, x_vaccini_tot(:,[1 13 25])) %visione suscettibili
legend('S(t)','S1(t)','S2(t)')

figure(3)
plot(t, x_vaccini_tot (:, 1)) %visione suscettibili 1 gatto

figure(4)
plot(t, x_vaccini_tot(:,2:12)) %altre variabili 1 gatto
legend('E1(t)','E2(t)', 'P(t)', 'I(t)', 'A(t)', 'H1(t)','H2(t)','H3(t)', 'Q(t)', 'R(t)', 'D(t)')

figure(5) %altre variabili 2 gatto
plot(t, x_vaccini_tot(:,14:24))
legend('E1(t)','E2(t)', 'P(t)', 'I(t)', 'A(t)', 'H1(t)','H2(t)','H3(t)', 'Q(t)', 'R(t)', 'D(t)')

figure(6) %altre variabili 3 gatto
plot(t, x_vaccini_tot(:,26:31))
legend('E21(t)','E22(t)', 'P2(t)', 'I2(t)', 'A2(t)','R2(t)')

%% rifacciamo l'analisi del tasso di raddoppio nel nuovo modello per vedere se c'� stato un rallentamento 
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
t_radd1 = t(selez1)

%qui addirittura sembra essere rilevante solo se consideriamo tutti i
%reparti


%% ora facciamo la cascata di 3 ode su infetti come suggerito da manfredi
options = odeset('RelTol',1e-8,'AbsTol',1e-10);
global x0_casc_inf
Lvect= [zeros(nolockdown,1); 0.6*ones(novax,1); 0.3*ones(136,1)]; %si suppone che il lockdown
%duri a partire dal 18esimo giorno e termini il giorno in cui cominciano le
%vaccinazioni
x0_casc_inf = [1-1*E0 , 1*E0, zeros(1,28)];
tic
[t,x_vaccini_tot]= ode45('gatto_vaccini_unico_cascatesoloInfetti', time, x0_casc_inf,options); 
toc

figure(1)
plot(t, x_vaccini_tot(:,:)) %plot totale

figure(2)
plot(t, x_vaccini_tot(:,[1 12 23])) %visione suscettibili
legend('S(t)','S1(t)','S2(t)')

figure(3)
plot(t, x_vaccini_tot(:,[2,3, 6:8])) %altre variabili 1 gatto
legend('E(t)','P(t)','I3(t)', 'A(t)', 'H(t)', 'Q(t)', 'R(t)', 'D(t)')

figure(4) %altre variabili 2 gatto
plot(t, x_vaccini_tot(:,13:22))
legend('E1(t)','P1(t)', 'I13(t)', 'A1(t)', 'H1(t)', 'Q1(t)', 'R1(t)', 'D1(t)')

figure(5) %altre variabili 3 gatto
plot(t, x_vaccini_tot(:,23:30))
legend('E2(t)','P2(t)', 'I12(t)','I22(t)','I23(t)', 'A2(t)', 'R2(t)')

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

