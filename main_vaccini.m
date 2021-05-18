clc
clear
clear all
%%
format longg
parameters_vaccini;
dati_vaccini;
global lambda deltaE deltaP sigma eta gammaI alfaI gammaA zeta gammaH ...
       alfaH gammaQ gammaA x0 N eff1 eff2 ef1 prima_dose_ seconda_dose_

% prima_d=prima_dose_norm;
% seconda_d=seconda_dose_norm;

%% Caricamento dati iniziali opzione 1
%si fa evolvere l'epidemia senza alcun controllo (lockdown)
%e senza usare vaccini per un numero specificato di giorni, poi si comincia
%l'intervento con i vaccini con le somministrazioni giornaliere effettuate
%in Italia a partire dal 27/12/2020
novax=80;
N=novax+136;
time= 0:1:N-1;
E0= 10/pop;
x0= [1-E0 E0 0 0 0 0 0 0 0]; 
prima_dose_= [zeros(novax,1); prima_dose_norm];
seconda_dose_=[zeros(novax,1); seconda_dose_norm];


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
%% modello Gatto con vaccini
tic
[t,x_vaccini]= ode45('gatto_vaccini', time, x0); 
toc

figure(1)
plot(t, x_vaccini(:,:))
legend('S(t)', 'E(t)', 'P(t)', 'I(t)', 'A(t)', 'H(t)', 'Q(t)', 'R(t)', 'D(t)')

figure(2)
plot(t, x_vaccini(:,2:end))
legend('E(t)', 'P(t)', 'I(t)', 'A(t)', 'H(t)', 'Q(t)', 'R(t)', 'D(t)')

% x10=x_vaccini(novax,:); % ho provato a mettere le condizioni iniziali per il secondo set di equazioni ... 
                        % come i valori della dinamica nel giorno in cui si
                        % comincia la vaccinazione
%% modello Gatto 2 popolazione (prima dose)
x10 = 1e-9*ones(9,1);
tic
[t,x_vaccini2]= ode45('gatto_vaccini2', time, x10); 
toc
%
figure(3)
plot(t, x_vaccini2(:,:))
legend('S(t)','E(t)', 'P(t)', 'I(t)', 'A(t)', 'H(t)', 'Q(t)', 'R(t)', 'D(t)')
figure(4)
plot(t, x_vaccini2(:,2:end))
legend('E(t)', 'P(t)', 'I(t)', 'A(t)', 'H(t)', 'Q(t)', 'R(t)', 'D(t)')
%% modello Gatto 3 popolazione (seconda dose)
x20 = 1e-9*ones(6,1);
tic
[t,x_vaccini2]= ode45('gatto_vaccini3', time, x20); 
toc
%
figure(5)
plot(t, x_vaccini2(:,:))
legend('S(t)','E(t)', 'P(t)', 'I(t)', 'A(t)','R(t)')

figure(6)
plot(t, x_vaccini2(:,2:end))
legend('E(t)', 'P(t)', 'I(t)', 'A(t)', 'R(t)')