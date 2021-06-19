clc
clear
clear all
%%
global lambda deltaE deltaP sig eta gammaI alfaI gammaA zeta gammaH ...
       alfaH gammaQ gammaA x0 N eff1 eff2 ef1 prima_d seconda_d Lvect NV R0
R0 = 3.6;
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
x0= [1-1*E0; 1*E0; zeros(22,1)];
tic
[x_vaccini_tot]= ode4(@gatto_vaccini_unico, 0, 1, nolockdown+novax-1, x0); 
toc
% lambda in questo file accoppia le sottodinamiche dei 3 set di equazioni dei vaccini
soluzione= zeros(novax+nolockdown, 24);
for j=1:1:295
    for i= 1:1:24
        soluzione(j,i)=x_vaccini_tot(24*(j-1)+i);
    end
end

figure(1)
plot(soluzione(:,:))%plot totale
legend('S(t)','E(t)','P(t)','I(t)','A(t)','H(t)','Q(t)','R(t)','D(t)')
xlabel('Days')
%lavorando su nolockdown+novax giorni non subentrano i vaccini
%% Cambio di R0
%si fa evolvere l'epidemia senza controllo per novax+nolockdown giorni in corrispondenza di nolockdown
%si inserisce il cambio di R0, poi si comincial'intervento con i vaccini con le somministrazioni
%giornaliere effettuate in Italia a partire dal 27/12/2020
% fonte: https://github.com/pcm-dpc/COVID-19/blob/master/dati-andamento-nazionale/dpc-covid19-ita-andamento-nazionale-20201227.csv
% fase iniziale (nolockdown giorni) con R0=3.6
Lvect = zeros(N,1);
x0= [1-1*E0; 1*E0; zeros(22,1)];
tic
[x_vaccini_tot]= ode4(@gatto_vaccini_unico, 0,1,nolockdown-1, x0); 
toc
%
soluzione1= zeros(nolockdown, 24);
for j=1:1:nolockdown
    for i= 1:1:24
        soluzione1(j,i)=x_vaccini_tot(24*(j-1)+i);
    end
end


figure(1)
plot(soluzione1(:,1))%plot totale
legend('S(t)','E(t)','P(t)','I(t)','A(t)','H(t)','Q(t)','R(t)','D(t)')
xlabel('Days')

figure(2)
plot( soluzione1(:,2:9))%plot totale
legend('E(t)','P(t)','I(t)','A(t)','H(t)','Q(t)','R(t)','D(t)')
xlabel('Days')

R0= 1.0; % si può settare l'R0 desiderato, r0_raddoppio aggiorna i beta_i
r0_raddoppio; %ricalcolo dei beta_i con R0=1.1
% carichiamo i nuovi parametri dei vaccini (si potrebbe fare anche con il lockdown, ma così scaliamo direttamente R0 e i beta)

x01 = soluzione1(nolockdown,:); %aggiorniamo le condizioni iniziali del nuovo sistema

tic
[x_vaccini_tot]= ode4(@gatto_vaccini_unico, 0,1,novax-1, x01'); 
toc
%
soluzione2= zeros(novax, 24);
for j=1:1:novax
    for i= 1:1:24
        soluzione2(j,i)=x_vaccini_tot(24*(j-1)+i);
    end
end

figure(3)
plot( soluzione2(:,1))%plot suscettibili
legend('S(t)')
xlabel('Days')

figure(4)
plot( soluzione2(:,2:9))%plot totale
legend('E(t)','P(t)','I(t)','A(t)','H(t)','Q(t)','R(t)','D(t)')
xlabel('Days')

figure(5)
soluzione= [soluzione1(:,:); soluzione2(:,:)];
plot(soluzione(:,:))
legend('S(t)','E(t)','P(t)','I(t)','A(t)','H(t)','Q(t)','R(t)','D(t)')
xlabel('Days')

%usiamo questi dati come condizione iniziale per l'ottimizzazione
x0_opt=soluzione(novax+nolockdown,:);
%% Introduzione del lockdown
%lockdown introdotto dopo nolockdown giorni
%questo caso equivale alla scalatura di R0 ma viene effettuata attraverso
%il lockdown: i beta_i non cambiano direttamente ma sono limitati dal fatto
%che settiamo un coefficiente (1-teta*L)^2 di fronte a lambda.
%e quindi riassegnamo a R0 il valore originale
R0 = 3.6;
Lvect = [zeros(nolockdown,1);0.635*ones(N-nolockdown,1)];
parameters_vaccini;
x0= [1-1*E0; 1*E0; zeros(22,1)];

tic
[x_vaccini_tot]= ode4(@gatto_vaccini_unico, 0, 1, novax+nolockdown-1, x0); 
toc
% lambda in questo file accoppia le sottodinamiche dei 3 set di equazioni dei vaccini
soluzione= zeros(novax+nolockdown, 24);
for j=1:1:novax+nolockdown
    for i= 1:1:24
        soluzione(j,i)=x_vaccini_tot(24*(j-1)+i);
    end
end

figure(1)
plot(soluzione(:,:))%plot totale
legend('S(t)','E(t)','P(t)','I(t)','A(t)','H(t)','Q(t)','R(t)','D(t)')
xlabel('Days')

figure(2)
plot(soluzione(:,2:9)) %primo set di equazioni
legend('E(t)', 'P(t)', 'I(t)', 'A(t)', 'H(t)', 'Q(t)', 'R(t)', 'D(t)')
xlabel('Days')

figure(3)
plot(soluzione(:,1))
legend('S(t)')
xlabel('Days')
%% ODE waterfall
%ora introduciamo le cascate di Ode ma secondo il Gatto (su E e su H come a pagina 11)
global x0_casc

parameters_vaccini;

Lvect= [zeros(nolockdown,1); 0*ones(novax,1)]; %da qua possiamo gestire il lockdown sui vari intervalli
x0_casc = [1-1*E0; 1*E0; zeros(29,1)]; %aumentano le CI aumentanto i compartimenti

tic
[x_vaccini_tot]= ode4(@gatto_vaccini_unico_cascate, 0, 1, nolockdown+novax-1, x0_casc); 
toc
% lambda in questo file accoppia le sottodinamiche dei 3 set di equazioni dei vaccini
soluzione= zeros(novax+nolockdown, 31);
for j=1:1:295
    for i= 1:1:31
        soluzione(j,i)=x_vaccini_tot(31*(j-1)+i);
    end
end

figure(1)
plot( soluzione(:,:))%plot totale
legend('S(t)','E(t)','P(t)','I(t)','A(t)','H(t)','Q(t)','R(t)','D(t)')
xlabel('Days')

figure(2)
plot(soluzione(:,[1 13 25])) %visione suscettibili (si apprezzano chiaramente solo quelli del primo set)
legend('S(t)','S1(t)','S2(t)')
xlabel('Days')

figure(3)
plot(soluzione(:,2:12)) %altre variabili 1 gatto
legend('E1(t)','E2(t)', 'P(t)', 'I(t)', 'A(t)', 'H1(t)','H2(t)','H3(t)', 'Q(t)', 'R(t)', 'D(t)')
xlabel('Days')

%% ODE waterfall 2
% ora facciamo la cascata di 3 ode su infetti come suggerito da manfredi
% (su infetti, esposti, asintomatici)

global x0_casc_inf

Lvect= [zeros(nolockdown,1); 0*ones(novax,1)]; %da qua possiamo gestire il lockdown sui vari intervalli
x0_casc_inf = [1-1*E0; 1*E0; zeros(40,1)];

tic
[x_vaccini_tot]= ode4(@gatto_vaccini_unico_cascatesoloInfetti, 0, 1, nolockdown+novax-1, x0_casc_inf); 
toc
% lambda in questo file accoppia le sottodinamiche dei 3 set di equazioni dei vaccini
soluzione= zeros(novax+nolockdown, 42);
for j=1:1:295
    for i= 1:1:42
        soluzione(j,i)=x_vaccini_tot(42*(j-1)+i);
    end
end

figure(1)
plot( soluzione(:,:))%plot totale
legend('S(t)','E(t)','P(t)','I(t)','A(t)','H(t)','Q(t)','R(t)','D(t)')
xlabel('Days')

figure(2)
plot(soluzione(:,1)) %visione suscettibili
legend('S(t)','S1(t)','S2(t)')

figure(3)
plot(soluzione(:,[2,3, 6:8])) %altre variabili 1 gatto
legend('E(t)','P(t)','I3(t)', 'A(t)', 'H(t)', 'Q(t)', 'R(t)', 'D(t)')

x0_casc_inf_opt=soluzione(novax+nolockdown,:);

%% Prova del funzionale di costo

global N_ott r ts xi w t_ott f R0
N_ott = N-nolockdown-novax; %lavoro su 165, con lockdown ottimale e vaccini
prima_d = prima_dose_norm;
seconda_d = seconda_dose_norm;

r=0.05;
xi=0; % termine aggiuntivo come extra costo delle vite
w=65000;
t_ott = 0:1:N_ott-1;
f= 0.99; % parametro per gestire lo sbilanciamento del funzionale (costo delle vite vs perdite economiche)
%siccome stiamo parlando del rilascio di lockdown, non ci interessa partire
%con valore di L a 0
x0= x0_opt; % condizioni iniziali per risolvere le ode
Lvect = zeros(N_ott,1);
U0= [0.65*ones((N_ott+1)/2,1);0.45*ones((N_ott+1)/2-1,1)]; %condizioni iniziali di vettore di ingresso di lockdown per l'ottimizzatore

lb= [0*ones((N_ott+1)/2,1);0*ones((N_ott+1)/2-1,1)]; % lower bounds
ub= 0.95*ones(N_ott,1); % upper bounds

%aggiorniamo il modello affinchè l'epidemia viaggi con tempo di raddoppio
%circa 3, settando qui R0=2.65
R0= 2.65;
r0_raddoppio; 

tic
[Uvec,fval,exitflag] = fmincon(@cost_function,U0,[],[],[],[],lb,ub,@nonlincon, options_lockdown);
toc
t = 0:N_ott-1;
% evolution 
Lvect = Uvec;
tic
[x_vaccini_tot]= ode4(@gatto_vaccini_unico, 0, 1, N_ott-1, x0_opt'); 
toc

soluzione= zeros(N_ott, 24);
for j=1:1:N_ott
    for i= 1:1:24
        soluzione(j,i)=x_vaccini_tot(24*(j-1)+i);
    end
end

figure(1)
plot(Uvec)
legend('Optimizated lockdown')
xlabel('Days')

figure(2)
plot(soluzione(:,:))
ylabel('22 variables')
xlabel('Days')

figure(3)
plot(soluzione(:, [2 3 4 5 6 7 8 9]))
xlabel('Days')
ylabel('No vaxed population')

%% %% Ottimizzazione con cascate di ODE
% 
% global N_ott r ts xi w t_ott f R0
% N_ott = N-nolockdown-novax; %lavoro su 165, con lockdown ottimale e vaccini
% prima_d = prima_dose_norm;
% seconda_d = seconda_dose_norm;
% 
% r=0.05;
% xi=0; % termine aggiuntivo come extra costo delle vite
% w=65000;
% t_ott = 0:1:N_ott-1;
% f= 0; % parametro per gestire lo sbilanciamento del funzionale (costo delle vite vs perdite economiche)
% %siccome stiamo parlando del rilascio di lockdown, non ci interessa partire
% %con valore di L a 0
% x0= x0_opt; % condizioni iniziali per risolvere le ode
% Lvect = zeros(N_ott,1);
% U0 = [0.5*ones(N_ott,1)]; %condizioni iniziali di vettore di ingresso di lockdown per l'ottimizzatore
% 
% lb= [0*ones((N_ott+1)/2,1);0*ones((N_ott+1)/2-1,1)]; % lower bounds
% ub= 0.95*ones(N_ott,1); % upper bounds
% 
% %aggiorniamo il modello affinchè l'epidemia viaggi con tempo di raddoppio
% %circa 3, settando qui R0=2.65
% R0= 2.65;
% r0_raddoppio; 
% 
% tic
% [Uvec,fval,exitflag] = fmincon(@cost_function_waterfall,U0,[],[],[],[],lb,ub,@nonlincon_casc, options_lockdown);
% toc
% t = 0:N_ott-1;
% % evolution 
% Lvect = Uvec;
% tic
% [x_vaccini_tot]= ode4(@gatto_vaccini_unico_cascatesoloInfetti, 0, 1, N_ott-1, x0'); 
% toc
% 
% soluzione= zeros(N_ott, 42);
% for j=1:1:N_ott
%     for i= 1:1:42
%         soluzione(j,i)=x_vaccini_tot(42*(j-1)+i);
%     end
% end
% 
% figure(1)
% plot(Uvec)
% legend('Optimizated lockdown')
% xlabel('Days')
% 
% figure(2)
% plot(soluzione(:,:))
% ylabel('42 variables')
% xlabel('Days')
% 
% figure(3)
% plot(soluzione(:, [2 3 4 5 6 7 8 9 10 11 12 13 14 15])) %primo set di equazioni


%

%% ODE waterfall con cambio di R0 
% ora facciamo la cascata di 3 ode su infetti come suggerito da manfredi
% (su infetti, esposti, asintomatici), per fornire le condizioni iniziale
% dopo 25 giorni nolockdown e 270 novax con R0=1.0

%cambiamo i parametri equivalenti di R0, poichè l'epidemia rallenta
%a causa delle nuove distribuzioni gamma
parameters_vaccini
dati_vaccini
prima_d= [zeros(novax+nolockdown,1); prima_dose_norm];
seconda_d=[zeros(novax+nolockdown,1); seconda_dose_norm];

global x0_casc_inf
R0=5.5;
r0_raddoppio;
Lvect= [zeros(nolockdown,1); zeros(novax,1)]; %da qua possiamo gestire il lockdown sui vari intervalli
x0_casc_inf = [1-1*E0; 1*E0; zeros(40,1)];

tic
[x_vaccini_tot]= ode4(@gatto_vaccini_unico_cascatesoloInfetti, 0, 1, nolockdown-1, x0_casc_inf); 
toc
% lambda in questo file accoppia le sottodinamiche dei 3 set di equazioni dei vaccini
soluzione1= zeros(nolockdown, 42);
for j=1:1:nolockdown
    for i= 1:1:42
        soluzione1(j,i)=x_vaccini_tot(42*(j-1)+i);
    end
end

R0= 1.1; % si può settare l'R0 desiderato, r0_raddoppio aggiorna i beta_i
r0_raddoppio; %ricalcolo dei beta_i con R0=1.0
x01_casc_inf= soluzione1(nolockdown,:)';

tic
[x_vaccini_tot]= ode4(@gatto_vaccini_unico_cascatesoloInfetti, 0, 1, novax-1, x01_casc_inf); 
toc
% lambda in questo file accoppia le sottodinamiche dei 3 set di equazioni dei vaccini
soluzione2= zeros(novax, 42);
for j=1:1:novax
    for i= 1:1:42
        soluzione2(j,i)=x_vaccini_tot(42*(j-1)+i);
    end
end

soluzione=[soluzione1(:,:);soluzione2(:,:)];

figure(1)
plot( soluzione(:,1))%plot totale
legend('S(t)','E(t)','P(t)','I(t)','A(t)','H(t)','Q(t)','R(t)','D(t)')
xlabel('Days')

x0_casc_inf_opt = soluzione(novax+nolockdown,:)';

%% Ottimizzazione con cascate di ODE

global N_ott r ts xi w t_ott f R0
N_ott = N-nolockdown-novax; %lavoro su 165, con lockdown ottimale e vaccini
prima_d = prima_dose_norm;
seconda_d = seconda_dose_norm;

r=0.05;
xi=0; % termine aggiuntivo come extra costo delle vite
w=65000;
t_ott = 0:1:N_ott-1;
f= 0.9; % parametro per gestire lo sbilanciamento del funzionale (costo delle vite vs perdite economiche)
%siccome stiamo parlando del rilascio di lockdown, non ci interessa partire
%con valore di L a 0
x0_casc_inf= x0_casc_inf_opt; % condizioni iniziali per risolvere le ode
Lvect = zeros(N_ott,1);
U0 = [0.7*ones(N_ott,1)]; %condizioni iniziali di vettore di ingresso di lockdown per l'ottimizzatore

lb= [0*ones((N_ott+1)/2,1);0*ones((N_ott+1)/2-1,1)]; % lower bounds
ub= 0.95*ones(N_ott,1); % upper bounds
%
%aggiorniamo il modello affinchè l'epidemia viaggi con tempo di raddoppio
%circa 3, settando qui R0=2.65
R0= 2.65;
r0_raddoppio; 

tic
[Uvec,fval,exitflag] = fmincon(@cost_function_waterfall,U0,[],[],[],[],lb,ub,@nonlincon_casc, options_lockdown);
toc
t = 0:N_ott-1;
% evolution 
Lvect = Uvec;
tic
[x_vaccini_tot]= ode4(@gatto_vaccini_unico_cascatesoloInfetti, 0, 1, N_ott-1, x0_casc_inf_opt); 
toc

soluzione= zeros(N_ott, 42);
for j=1:1:N_ott
    for i= 1:1:42
        soluzione(j,i)=x_vaccini_tot(42*(j-1)+i);
    end
end

figure(8)
plot(Uvec)
legend('Optimizated lockdown')
xlabel('Days')

figure(9)
plot(soluzione(:,:))
ylabel('42 variables')
xlabel('Days')


%%
% %% ora proviamo ad impostare l'ottimizzazione parametrica: tipo lockdown costante su tutta la finestra temporale
% %come prima aggiustando R0 per avere il tasso di raddoppio richiesto
% 
% %anche qui scalare nellos script R0 ai valori corrispondenti al tasso di
% %raddoppio voluto
% R0=2,65;
% r0_raddoppio;
% 
% global N_ott  r ts xi w t_ott
% N_ott = N-nolockdown-novax; %che è 165
% prima_d = prima_dose_norm;
% seconda_d = seconda_dose_norm;
% 
% r=0.05;
% xi=0;
% w=65000;
% t_ott = 0:1:N_ott-1;
% 
% valori0 = 0.9;
% 
% %costruisco un vettore di costanti tale da avere lockdownccostanti su 14 giorni
% for i = 1:14:N_ott
%    valori0 = [valori0 , 0.9]; 
% end
% 
% % %prova eventuale scalatura di valori0
% % finestre = 1:length(valori0);
% % scala = (length(valori0)- finestre)/length(valori0);
% % valori0 = valori0.*scala;
% 
% % plot(valori0)
% 
% %
% n = length(valori0);
% lb= -0.00001*ones(n,1); % lower bounds lo mettiamo negativo leggermente 
% %per fare in modo che sia verificata strettamente la disuguglianza
% 
% ub= 0.9.*ones(n,1); % upper bounds
% 
% tic
% [Uvec_param,fval,exitflag] = fmincon('cost_function_param',valori0,[],[],[],[],lb,ub,[]);
% toc
% %
% % evolution 
% Lvect = [];
% for i = 1:n;
%  Lvect = [Lvect,Uvec_param(i)*ones(1,14)]   ; 
% end
% %
% % t = 0:N_ott-1 ;
% [XFin] = ode4(@gatto_vaccini_unico,0,1,N_ott-1,x0);
% % 
% figure(1)
% plot(Lvect(1:N_ott))
% xlabel('Days')
% 
% 
% tic
% [x_vaccini_tot]= ode4(@gatto_vaccini_unico, 0, 1, N_ott-1, x0' ); 
% toc
% % lambda in questo file accoppia le sottodinamiche dei 3 set di equazioni dei vaccini
% soluzione= zeros(N_ott, 24);
% for j=1:1:165
%     for i= 1:1:24
%         soluzione(j,i)=x_vaccini_tot(24*(j-1)+i);
%     end
% end
% %%
% figure(2)
% plot( soluzione(:,1))%plot totale
% legend('S(t)','E(t)','P(t)','I(t)','A(t)','H(t)','Q(t)','R(t)','D(t)')
% xlabel('Days')
% 
% 
% 
% 
% %% logistiche 
% % optimal control
% 
% U0 = [0.6 0.8 30 ... 
%     0.6 0.4 100]; 
% lb = [0 0 0 0 0 0]; 
% ub = [0.9 1 N-1 0.9 1 N-1];
% 
% %lb= zeros(N,1); % lower bounds
% %ub= 0.9.*ones(N,1); % upper bounds
% 
% 
% options = optimoptions('fmincon','Display','iter-detailed');
% tic
% [Uvec,fval,exitflag] =fmincon('cost_function_param_logi',U0,[],[],[],[],lb,ub,'nonlincon',options);
% toc
% t = 0:N_ott-1;
% % evolution 
% %Lvect= [zeros(nolockdown,1); 0.6*ones(novax,1); 0.3*ones(NV,1)];
% Lvect = Utime2par(Uvec,t);
% %Lvect= Uvec;
% [time, XFin] = ode45('gatto_vaccini_unico',t,x0);
% 
% %plots
% figure
% t5 = tiledlayout(5,1);
% nexttile; plot(time,Lvect.*100,'r'); ylabel('Lockdown (%)'); ylim([0 100]);
% nexttile; plot(time,XFin(:,1).*100); hold on; plot(time,XFin(:,10).*100,'r');hold on; plot(time,XFin(:,19).*100,'b');
% ylabel('S (%)')
% nexttile; plot(time,XFin(:,2).*100); hold on; plot(time,XFin(:,11).*100,'r');hold on; plot(time,XFin(:,20).*100,'b');
% ylabel('E (%)')
% nexttile; plot(time,XFin(:,8).*100); hold on; plot(time,XFin(:,17).*100,'r');hold on; plot(time,XFin(:,24).*100,'b');
% ylabel('R (%)')
% nexttile; plot(time,100-XFin(:,9).*100);  hold on; plot(time,100-XFin(:,18).*100,'r');
% ylabel('D (%)')
% title(t5,'Optimal Lockdown'); xlabel(t5,'Time (days)');
% legend({'No lockdown' 'Optimal lockdown'},'orientation','horizontal','location','southoutside');

