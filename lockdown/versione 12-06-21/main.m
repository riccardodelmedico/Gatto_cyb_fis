clc
clear
clear all
%%
global lambda deltaE deltaP sig eta gammaI alfaI gammaA zeta gammaH ...
       alfaH gammaQ gammaA x0 N eff1 eff2 ef1 prima_d seconda_d Lvect NV
   
format longg %utilizzando valori molto piccoli le 4 cifre decimali non sono sufficienti
dati_vaccini;
parameters_vaccini;

time= 0:1:N-1; % N viene definitio in parameters_vaccini
E0= 10/pop; % condizione iniziale sul numero di esposti
 
prima_d= [zeros(novax+nolockdown,1); prima_dose_norm];
seconda_d=[zeros(novax+nolockdown,1); seconda_dose_norm];
% si suppone di avere novax giorni senza vaccini
figure(1)
plot(time, prima_d)
hold on
plot(time,seconda_d)
hold off

legend('First dose','Second dose')
xlabel('Days')
ylabel('Vaccines per days')
%% Dinamiche accoppiate, senza vaccinazione (Gatto model)
%si fa evolvere l'epidemia senza controllo per novax giorni, poi si comincia
%l'intervento con i vaccini con le somministrazioni giornaliere effettuate
%in Italia a partire dal 27/12/2020
% fonte: https://github.com/pcm-dpc/COVID-19/blob/master/dati-andamento-nazionale/dpc-covid19-ita-andamento-nazionale-20201227.csv

Lvect= zeros(N,1); % prima soluzione in evoluzione libera
x0= [1-1*E0 1*E0 zeros(1,22)];

tic
[t0,x_vaccini_tot]= ode45('gatto_vaccini_unico', 0:1:nolockdown+novax-1, x0,options_ode); 
toc
% lambda in questo file accoppia le sottodinamiche dei 3 set di equazioni dei vaccini
figure(1)
plot(t0, x_vaccini_tot(:,:))%plot totale
legend('S(t)','E(t)','P(t)','I(t)','A(t)','H(t)','Q(t)','R(t)','D(t)')
xlabel('Days')
% lavorando su nolockdown+novax giorni non subentrano i vaccini
%% Cambio di R0
%si fa evolvere l'epidemia senza controllo per novax+nolockdown giorni in corrispondenza di nolockdown
%si inserisce il cambio di R0, poi si comincial'intervento con i vaccini con le somministrazioni
%giornaliere effettuate in Italia a partire dal 27/12/2020
% fonte: https://github.com/pcm-dpc/COVID-19/blob/master/dati-andamento-nazionale/dpc-covid19-ita-andamento-nazionale-20201227.csv
% fase iniziale (nolockdown giorni) con R0=3.6
Lvect = zeros(N,1);
x0= [1-1*E0 1*E0 zeros(1,22)];
tic
[~,x_vaccini_tot]= ode45('gatto_vaccini_unico', 0:1:nolockdown-1, x0,options_ode); 
toc

parameters_vaccini_R0_1; %ricalcolo dei beta_i con R0=1.1
% carichiamo i nuovi parametri dei vaccini (si potrebbe fare anche con il lockdown, ma così scaliamo direttamente R0 e i beta)

x01 = x_vaccini_tot(end,:); %aggiorniamo le condizioni iniziali del nuovo sistema

%risolviamo il Gatto senza vaccini per 9 mesi (novax giorni), poi intervengono i vaccini
tic
[~,x_vaccini_tot1]= ode45('gatto_vaccini_unico', 0:1:novax+NV-1, x01,options_ode); 
toc

x_vaccini=[x_vaccini_tot(:,:); x_vaccini_tot1(:,:)]; %uniamo i risultati delle due soluzioni
figure(1)
plot( [x_vaccini(:,[1 10 19])]) %visione suscettibili
legend('S(t)','S1(t)','S2(t)')
xlabel('Days')
%
figure(2)
plot(x_vaccini(:,2:9)) %primo set di equazioni
legend('E(t)', 'P(t)', 'I(t)', 'A(t)', 'H(t)', 'Q(t)', 'R(t)', 'D(t)')
xlabel('Days')
%
figure(3) %altre variabili: secondo set di equazioni
plot(x_vaccini(:,11:18))
legend('E1(t)', 'P1(t)', 'I1(t)', 'A1(t)', 'H1(t)', 'Q1(t)', 'R1(t)', 'D1(t)')
xlabel('Days')
% 
figure(4) %altre variabili 3 gatto
plot(x_vaccini(:,20:24))
legend('S2(t)','E2(t)', 'P2(t)', 'I2(t)', 'A2(t)','R2(t)')
xlabel('Days')

%% Introduzione del lockdown
%lockdown introdotto dopo nolockdown giorni
%questo caso equivale alla scalatura di R0 ma viene effettuata attraverso
%il lockdown: i beta_i non cambiano direttamente ma sono limitati dal fatto
%che settiamo un coefficiente (1-teta*L)^2 di fronte a lambda.
Lvect = [zeros(nolockdown,1);0.8*ones(N-nolockdown,1)];
parameters_vaccini;
x0= [1-1*E0 1*E0 zeros(1,22)];

tic
[~,x_vaccini_lockdown]= ode45('gatto_vaccini_unico',0:1:N-1,x0,options_ode);
toc

figure(1)
plot(x_vaccini_lockdown(:,1:9))
legend('S(t)','E(t)','P(t)','I(t)','A(t)','H(t)','Q(t)','R(t)','D(t)')
xlabel('Days')

figure(2)
plot(x_vaccini_lockdown(:,2:9)) %primo set di equazioni
legend('E(t)', 'P(t)', 'I(t)', 'A(t)', 'H(t)', 'Q(t)', 'R(t)', 'D(t)')
xlabel('Days')
%
figure(3) %altre variabili: secondo set di equazioni
plot(x_vaccini_lockdown(:,11:18))
legend('E1(t)', 'P1(t)', 'I1(t)', 'A1(t)', 'H1(t)', 'Q1(t)', 'R1(t)', 'D1(t)')
xlabel('Days')
% 
figure(4) %altre variabili 3 gatto
plot(x_vaccini_lockdown(:,20:24))
legend('S2(t)','E2(t)', 'P2(t)', 'I2(t)', 'A2(t)','R2(t)')
xlabel('Days')

%% ODE waterfall
%ora introduciamo le cascate di Ode ma secondo il Gatto (su E e su H come a pagina 11)
global x0_casc

parameters_vaccini;
dati_vaccini;

Lvect= [zeros(N,1)]; %lockdown nullo
x0_casc = [1-1*E0 , 1*E0, zeros(1,29)]; %aumentano le CI aumentanto i compartimenti

options = odeset('RelTol',1e-8,'AbsTol',1e-10);
tic
[t,x_vaccini_waterfall]= ode45('gatto_vaccini_unico_cascate', 0:1:novax+nolockdown-1, x0_casc,options); 
toc

figure(1)
plot(t, x_vaccini_waterfall(:,:)) %plot totale
xlabel('Days')

figure(2)
plot(t, x_vaccini_waterfall(:,[1 13 25])) %visione suscettibili (si apprezzano chiaramente solo quelli del primo set)
legend('S(t)','S1(t)','S2(t)')
xlabel('Days')

figure(3)
plot(t, x_vaccini_waterfall(:,2:12)) %altre variabili 1 gatto
legend('E1(t)','E2(t)', 'P(t)', 'I(t)', 'A(t)', 'H1(t)','H2(t)','H3(t)', 'Q(t)', 'R(t)', 'D(t)')
xlabel('Days')

%% ODE waterfall 2
% ora facciamo la cascata di 3 ode su infetti come suggerito da manfredi
% (su infetti, esposti, asintomatici)
options = odeset('RelTol',1e-8,'AbsTol',1e-10);
global x0_casc_inf
Lvect= [zeros(nolockdown,1); 0*ones(novax,1); 0*ones(NV,1)]; %da qua possiamo gestire il lockdown sui vari intervalli
x0_casc_inf = [1-1*E0 , 1*E0, zeros(1,40)];
time = 0:1:novax+nolockdown-1;
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


%% %% Prova del funzionale di costo

global N_ott;
N_ott = N-nolockdown-novax; %che è 165
prima_d = prima_dose_norm;
seconda_d = seconda_dose_norm;

global r ts xi w t_ott
r=0.05;
ts=1;
xi=0;
w=65000;
t_ott = 0:1:N_ott-1;

%siccome stiamo parlando del rilascio di lockdown, non ci interessa partire
%con valore di L a 0

Lvect = zeros(1,N_ott);
U0 = [0.7*ones(N_ott,1)];

lb= zeros(N_ott,1); % lower bounds
ub= 0.9.*ones(N_ott,1); % upper bounds

%ed attenzione che devo aggiornare il modello affinchè l'epidemia viaggi
%con tempo di raddoppio circa 3, settando qui R0 a 2.3
%
parameters_vaccini_R0_raddoppio;

tic
[Uvec,fval,exitflag] = fmincon('cost_function',U0,[],[],[],[],lb,ub,[],options_lockdown);
toc
t = 0:N_ott-1;
% evolution 
Lvect = Uvec;
[time, XFin] = ode45('gatto_vaccini_unico',t,x0);
%
figure(1)
plot(time, Uvec)

figure(2)
plot(time, XFin(:,:))

%% ora proviamo ad impostare l'ottimizzazione parametrica: tipo lockdown costante su tutta la finestra temporale
%come prima aggiustando R0 per avere il tasso di raddoppio richiesto

%anche qui scalare nellos script R0 ai valori corrispondenti al tasso di
%raddoppio voluto

parameters_vaccini_R0_raddoppio;

global N_ott;
N_ott = N-nolockdown-novax; %che è 165
prima_d = prima_dose_norm;
seconda_d = seconda_dose_norm;

global r ts xi w t_ott
r=0.05;
ts=1;
xi=0;
w=65000;
t_ott = 0:1:N_ott-1;

options_lockdown= optimoptions('fmincon','Display','iter-detailed','Algorithm','active-set',...
     'FunValCheck','on');

%time1=0:1:N-nolockdown-1;
% U0 = [zeros(nolockdown,1); 0.7.*ones(N-nolockdown,1).*(1-time1'./(N-nolockdown))];
%come sopra non avrebbe molto senso parlare di condizione con 0 lockdown
%per cui diamo come condizione iniziale un valore per cui Rt circa = 1
valori0 = 0.1;

%costruisco un vettore di costanti tale da avere lockdownccostanti su 14 giorni
for i = 1:14:N_ott
   valori0 = [valori0 , 0.1]; 
end

% %prova eventuale scalatura di valori0
% finestre = 1:length(valori0);
% scala = (length(valori0)- finestre)/length(valori0);
% valori0 = valori0.*scala;

% plot(valori0)

%
n = length(valori0);
lb= -0.00001*ones(n,1); % lower bounds lo mettiamo negativo leggermente 
%per fare in modo che sia verificata strettamente la disuguglianza

ub= 0.9.*ones(n,1); % upper bounds

tic
[Uvec_param,fval,exitflag] = fmincon('cost_function_param',valori0,[],[],[],[],lb,ub,[]);
toc
%
% evolution 
Lvect = [];
for i = 1:n;
 Lvect = [Lvect,Uvec_param(i)*ones(1,14)]   ; 
end
%
t = 0:N_ott-1 ;
[time, XFin] = ode45('gatto_vaccini_unico',t,x0);
% 
figure(1)
plot(time, Lvect(1:length(time)))

figure(2)
plot(time, XFin(:,:))


%% logistiche 
% optimal control

U0 = [0.6 0.8 30 ... 
    0.6 0.4 100]; 
lb = [0 0 0 0 0 0]; 
ub = [0.9 1 N-1 0.9 1 N-1];

%lb= zeros(N,1); % lower bounds
%ub= 0.9.*ones(N,1); % upper bounds


options = optimoptions('fmincon','Display','iter-detailed');
tic
[Uvec,fval,exitflag] =fmincon('cost_function_param_logi',U0,[],[],[],[],lb,ub,'nonlincon',options);
toc
t = 0:N_ott-1;
% evolution 
%Lvect= [zeros(nolockdown,1); 0.6*ones(novax,1); 0.3*ones(NV,1)];
Lvect = Utime2par(Uvec,t);
%Lvect= Uvec;
[time, XFin] = ode45('gatto_vaccini_unico',t,x0);

%plots
figure
t5 = tiledlayout(5,1);
nexttile; plot(time,Lvect.*100,'r'); ylabel('Lockdown (%)'); ylim([0 100]);
nexttile; plot(time,XFin(:,1).*100); hold on; plot(time,XFin(:,10).*100,'r');hold on; plot(time,XFin(:,19).*100,'b');
ylabel('S (%)')
nexttile; plot(time,XFin(:,2).*100); hold on; plot(time,XFin(:,11).*100,'r');hold on; plot(time,XFin(:,20).*100,'b');
ylabel('E (%)')
nexttile; plot(time,XFin(:,8).*100); hold on; plot(time,XFin(:,17).*100,'r');hold on; plot(time,XFin(:,24).*100,'b');
ylabel('R (%)')
nexttile; plot(time,100-XFin(:,9).*100);  hold on; plot(time,100-XFin(:,18).*100,'r');
ylabel('D (%)')
title(t5,'Optimal Lockdown'); xlabel(t5,'Time (days)');
legend({'No lockdown' 'Optimal lockdown'},'orientation','horizontal','location','southoutside');

