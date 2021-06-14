clc
clear
clear all
%%
global lambda deltaE deltaP sig eta gammaI alfaI gammaA zeta gammaH ...
       alfaH gammaQ gammaA x0 N eff1 eff2 ef1 prima_d seconda_d Lvect NV

   
format longg
dati_vaccini;
parameters_vaccini;



% prima_d=prima_dose_norm;
% seconda_d=seconda_dose_norm;

% Caricamento dati iniziali opzione 1
%si fa evolvere l'epidemia con controllo dopo 18 giorni(lockdown)
%e senza usare vaccini per un numero specificato di giorni, poi si comincia
%l'intervento con i vaccini con le somministrazioni giornaliere effettuate
%in Italia a partire dal 27/12/2020


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

dati_vaccini;
parameters_vaccini %i parametri sono quelli originali

options = odeset('RelTol',1e-7,'AbsTol',1e-8);
Lvect= zeros(N,1);% 0.5*ones(novax,1); 0.5*ones(NV,1)]; %si suppone che il lockdown
%duri a partire dal 25esimo giorno e termini il giorno in cui cominciano le
%vaccinazioni
x0= [1-1*E0 1*E0 zeros(1,22)];
%x0=[1-1*E0 1*E0 zeros(1,24)];

%risolviamo allora solo sull'intervallo di nolockdown, per cercare il tasso
%di raddoppio dell'epidemia libera

tic
[t,x_vaccini_tot]= ode45('gatto_vaccini_unico', 0:1:nolockdown+novax-1, x0,options); 
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

%% ora cambiamo R0, così scalando tutto e vediamo come evolve il sistema
% carichiamo i nuovi parametri dei vaccini 
%(si potrebbe fare anche con il lockdown, ma così scaliamo direttamente R0 e i beta)
tic
[t,x_vaccini_tot]= ode45('gatto_vaccini_unico', 0:1:nolockdown-1, x0,options); 
toc

figure(1)
plot(t, x_vaccini_tot(:,2:9)) %altre variabili 1 gatto
legend('E(t)', 'P(t)', 'I(t)', 'A(t)', 'H(t)', 'Q(t)', 'R(t)', 'D(t)')

parameters_vaccini_R0_1;
Lvect = zeros(N,1);

%aggiorniamo le condizioni iniziali del nuovo sistema
x01 = x_vaccini_tot(end,:);% variazione suscettibili

%risolviamo il Gatto senza vaccini per 6 mesi
tic
[t1,x_vaccini_tot1]= ode45('gatto_vaccini_unico', 0:1:novax, x01,options); 
toc

figure (2)
t_tot = 0:1:(novax+nolockdown);
plot(t_tot, [x_vaccini_tot(:,1); x_vaccini_tot1(:,1)]) %visione suscettibili

figure(3)
plot(t_tot, [x_vaccini_tot(:,2:end); x_vaccini_tot1(:,2:end)]) %visione restanti

%e aggiorniamo x0 per la soluzione del problema di controllo ottimo

x0 = x_vaccini_tot1(end, :);

 %% equivalente scalatura di R0, assunto che il lockdown abbia contributo quadratico
 %semplicemente vogliamo una scalatura di 3.6 del termine lambda, tramite
 %la moltiplicazione di (1-L)^2 = 1/3.6=0.2777---->(1-L)= 0.527--->L = 0.473
 
%rimettiamo allora R0 a 3.6
% parameters_vaccini;
% Lvect = [zeros(nolockdown,1); 0.48*ones(N-nolockdown,1) ];
% 
% %rimettiamo le condizioni iniziali come prima
% x0= [1-1*E0 1*E0 zeros(1,22)];
% %x0= [1-1*E0 1*E0 zeros(1,24)];
% %risolviamo il Gatto senza vaccini per tutti i mesi prima dei vaccini
% tic
% [t1,x_vaccini_tot2]= ode45('gatto_vaccini_unico', 0:1:(novax+nolockdown-1), x0,options); 
% toc
% 
% figure (1)
% t_tot = 0:1:(novax+nolockdown-1);
% plot(t_tot, x_vaccini_tot2(:,1))  % visione suscettibili
% 
% figure(2)
% plot(t_tot, x_vaccini_tot2(:,2:end))  % visione restanti variabili gatto


%% ora introduciamo le cascate di Ode ma secondo il Gatto (su E e su H come a pagina 11)
global x0_casc

parameters_vaccini;

Lvect= [zeros(nolockdown,1); 0.4*ones(novax,1); 0.4*ones(NV,1)]; %si suppone che il lockdown
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


%% ora facciamo la cascata di 3 ode su infetti come suggerito da manfredi
options = odeset('RelTol',1e-8,'AbsTol',1e-10);
global x0_casc_inf
Lvect= [zeros(nolockdown,1); 0.6*ones(novax,1); 0.3*ones(NV,1)]; %si suppone che il lockdown
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


%% %% Prova del funzionale di costo
%options_lockdown= optimoptions('fmincon','Display','none','Algorithm','active-set',...
 %   'OptimalityTolerance', 1e-5, 'MaxFunctionEvaluations', 5,'FunValCheck','on');
options_lockdown = optimoptions('fmincon','Display','iter-detailed','Algorithm','active-set','FunValCheck','on');
%time1=0:1:N-nolockdown-1;
% U0 = [zeros(nolockdown,1); 0.7.*ones(N-nolockdown,1).*(1-time1'./(N-nolockdown))];

global N_ott;
N_ott = N-nolockdown-novax; %che è 165
prima_d = prima_dose_norm;
seconda_d = seconda_dose_norm;


Lvect = zeros(1,N_ott);
U0 = [zeros(nolockdown,1); 0.7*ones(N_ott-nolockdown,1)];

lb= zeros(N_ott,1); % lower bounds
ub= 0.9.*ones(N_ott,1); % upper bounds
%
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

options_lockdown= optimoptions('fmincon','Display','iter-detailed','Algorithm','active-set',...
     'FunValCheck','on');

%time1=0:1:N-nolockdown-1;
% U0 = [zeros(nolockdown,1); 0.7.*ones(N-nolockdown,1).*(1-time1'./(N-nolockdown))];
valori0 = [0,0];

%costruisco un vettore di costanti tale da avere lockdownccostanti su 14 giorni
for i = 1:14:N_ott
   valori0 = [valori0 , 0.7]; 
end

%prova eventuale scalatura di valori0
finestre = 1:length(valori0);
scala = (length(valori0)- finestre)/length(valori0);
valori0 = valori0.*scala;

plot(valori0)


n = length(valori0);
lb= zeros(n,1); % lower bounds
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

