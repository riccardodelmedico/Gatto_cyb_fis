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
x0= [1-E0 E0 zeros(1,22)]; 
prima_d= [zeros(novax+nolockdown,1); prima_dose_norm];
seconda_d=[zeros(novax+nolockdown,1); seconda_dose_norm];
Lvect= [zeros(nolockdown,1); 0.7.*ones(novax,1); 0.45*ones(136,1)]; %si suppone che il lockdown
%duri a partire dal 18esimo giorno e termini il giorno in cui cominciano le
%vaccinazioni



% %% Caricamento dati iniziali opzione 2
% %si utilizzano i valori reali dell'epidemia nel giorno in cui sono
% %cominciate le vaccinazioni (27/12/2020) (dati da
% %https://github.com/pcm-dpc/COVID-19/blob/master/dati-andamento-nazionale/dpc-covid19-ita-andamento-nazionale-20201227.csv)
% %ipotesi: gli ospedalizzati comprendo TI e ricoverati
% %ipotesi: gli attualmente positivi sono divisi al 30% asintomatici,
% % 20% infettivi e 50% presintomatici
% N=136;
% time= 0:1:N-1;
% totale_ospedalizzati= 26151;
% totale_quarantena= 555609;
% totale_positivi= 581760;
% deceduti= 71925;
% dimessi_guariti= 1394011; %sarebbero i nostri recuperati








%%
%VEDIAMO CHE SUCCEDE ACCOPPIANDO LE DINAMICHE

tic
[t,x_vaccini_tot]= ode45('gatto_vaccini_unico', time, x0); 
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
plot(t, x_vaccini_tot(:,19:24))
legend('S2(t)','E2(t)', 'P2(t)', 'I2(t)', 'A2(t)','R2(t)')
%