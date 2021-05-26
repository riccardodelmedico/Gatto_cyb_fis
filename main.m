clc
clear
clear all
%%

format longg
parameters_vaccini;
dati_vaccini;
global lambda deltaE deltaP sig eta gammaI alfaI gammaA zeta gammaH ...
       alfaH gammaQ gammaA x0 N eff1 eff2 ef1 prima_d seconda_d Lvect NV

% prima_d=prima_dose_norm;
% seconda_d=seconda_dose_norm;

% Caricamento dati iniziali opzione 1
%si fa evolvere l'epidemia con controllo dopo 18 giorni(lockdown)
%e senza usare vaccini per un numero specificato di giorni, poi si comincia
%l'intervento con i vaccini con le somministrazioni giornaliere effettuate
%in Italia a partire dal 27/12/2020
nolockdown= 18;
novax= 62;
N= nolockdown + novax + NV;
time= 0:1:N-1;
E0= 100/pop;
 
prima_d= [zeros(novax+nolockdown,1); prima_dose_norm];
seconda_d=[zeros(novax+nolockdown,1); seconda_dose_norm];
%%
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
options = odeset('RelTol',1e-8,'AbsTol',1e-10);
Lvect= [zeros(nolockdown,1); 0.5*ones(novax,1); 0.5*ones(136,1)]; %si suppone che il lockdown
%duri a partire dal 18esimo giorno e termini il giorno in cui cominciano le
%vaccinazioni
x0= [1-1*E0 1*E0 zeros(1,22)];
tic
[t,x_vaccini_tot]= ode45('gatto_vaccini_unico', time, x0,options); 
toc

figure(1)
plot(t, x_vaccini_tot(:,:)) %plot totale

figure(2)
plot(t, x_vaccini_tot(:,[1 10 19])) %visione suscettibili
legend('S(t)','S1(t)','S2(t)')

figure(3)
plot(t, x_vaccini_tot(:,2:9)) %altre variabili 1 gatto
legend('E(t)', 'P(t)', 'I(t)', 'A(t)', 'H(t)', 'Q(t)', 'R(t)', 'D(t)')

figure(4) %altre variabili 2 gatto
plot(t, x_vaccini_tot(:,11:18))
legend('E1(t)', 'P1(t)', 'I1(t)', 'A1(t)', 'H1(t)', 'Q1(t)', 'R1(t)', 'D1(t)')

figure(5) %altre variabili 3 gatto
plot(t, x_vaccini_tot(:,20:24))
legend('S2(t)','E2(t)', 'P2(t)', 'I2(t)', 'A2(t)','R2(t)')
%

%% ora introduciamo le cascate di Ode ma secondo il Gatto (su E e su H come a pgina 11)
global x0_casc
Lvect= [zeros(nolockdown,1); 0.4*ones(novax,1); 0.3*ones(136,1)]; %si suppone che il lockdown
%duri a partire dal 18esimo giorno e termini il giorno in cui cominciano le
%vaccinazioni
x0_casc = [1-1*E0 , 1*E0, zeros(1,29)];
options = odeset('RelTol',1e-8,'AbsTol',1e-10);
tic
[t,x_vaccini_tot]= ode45('gatto_vaccini_unico_cascate', time, x0_casc,options); 
toc

figure(1)
plot(t, x_vaccini_tot(:,:)) %plot totale

figure(2)
plot(t, x_vaccini_tot(:,[1 13 25])) %visione suscettibili
legend('S(t)','S1(t)','S2(t)')

figure(3)
plot(t, x_vaccini_tot(:,2:12)) %altre variabili 1 gatto
legend('E1(t)','E2(t)', 'P(t)', 'I(t)', 'A(t)', 'H1(t)','H2(t)','H3(t)', 'Q(t)', 'R(t)', 'D(t)')

figure(4) %altre variabili 2 gatto
plot(t, x_vaccini_tot(:,14:24))
legend('E1(t)','E2(t)', 'P(t)', 'I(t)', 'A(t)', 'H1(t)','H2(t)','H3(t)', 'Q(t)', 'R(t)', 'D(t)')

figure(5) %altre variabili 3 gatto
plot(t, x_vaccini_tot(:,26:31))
legend('E21(t)','E22(t)', 'P2(t)', 'I2(t)', 'A2(t)','R2(t)')

%% ora facciamo la cascata di 3 ode su infetti come suggerito da manfredi
options = odeset('RelTol',1e-8,'AbsTol',1e-10);
global x0_casc_inf
Lvect= [zeros(nolockdown,1); 0.6*ones(novax,1); 0.5*ones(136,1)]; %si suppone che il lockdown
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

%% analisi andamento lambda
x = x_vaccini_tot;


S = x(:,1); %suscettibili
E = x(:,2); %esposti
P = x(:,3); %presintomatici
I1 = x(:,4);%infetti
I2 = x(:,5);
I3 = x(:,6);
A = x(:,7); %asintomatici
H = x(:,8); %ospedalizzati
Q = x(:,9); %quarantena casalinga
R = x(:,10); %recuperati
D = x(:,11); %dead
% parameters_vaccini;
% dati_vaccini;
S1 = x(:,12); %suscettibili
E1 = x(:,13); %esposti
P1 = x(:,14); %presintomatici
I11 = x(:,15);%infetti
I12 = x(:,16);
I13 = x(:,17);
A1 = x(:,18); %asintomatici
H1 = x(:,19); %ospedalizzati
Q1 = x(:,20); %quarantena casalinga
R1 = x(:,21); %recuperati
D1 = x(:,22); %dead
% parameters_vaccini;
% dati_vaccini;
S2 = x(:,23); %suscettibili
E2 = x(:,24); %esposti gruppo
P2 = x(:,25); %presintomatici
I21 = x(:,26);%infetti
I22 = x(:,27);
I23 = x(:,28);
A2 = x(:,29); %asintomatici
R2 = x(:,30); 

lambda= (betaP*(P+P1+P2) + betaI*(I3+I13+I23) + betaA*(A+A1+A2))./(S + E + (I3+I13+I23) + (A+A1+A2) + R);
plot(time,lambda)
%
plot(time,E)
