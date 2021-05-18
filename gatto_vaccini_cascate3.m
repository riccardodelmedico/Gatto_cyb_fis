clc
clear
clear all
%%
global lambda deltaE deltaP sigma eta gammaI alfaI gammaA zeta gammaH ...
       alfaH gammaQ gammaA Lvect x0 N eff1 eff2 ef1 prima_d seconda_d x10 x20
parameters_vaccini;
dati_vaccini;

prima_d=prima_dose_norm;
seconda_d=seconda_dose_norm;
format longg
N=136;
time= 0:1:N-1;
E0= 10/pop;
x0= [1-E0 E0 0 0 0 0 0 0 0]; 
%% modello Gatto con vaccini
tic
[t,x_vaccini]= ode45('gatto_vaccini', time, x0); 
toc
%
figure(1)
plot(t, x_vaccini(:,:))
legend('S(t)','E(t)', 'P(t)', 'I(t)', 'A(t)', 'H(t)', 'Q(t)', 'R(t)', 'D(t)')
figure(2)
plot(t, x_vaccini(:,2:end))
legend('E(t)', 'P(t)', 'I(t)', 'A(t)', 'H(t)', 'Q(t)', 'R(t)', 'D(t)')
% modello Gatto 2 popolazione (prima dose)
x10 = 1e-9*ones(9,1);
tic
[t,x_vaccini2]= ode45('gatto_vaccini2', time, x10); 
toc
%
figure(3)
plot(t, x_vaccini2(:,:))
legend('S1(t)','E1(t)', 'P1(t)', 'I1(t)', 'A1(t)', 'H1(t)', 'Q1(t)', 'R1(t)', 'D1(t)')
figure(4)
plot(t, x_vaccini2(:,2:end))
legend('E1(t)', 'P1(t)', 'I1(t)', 'A1(t)', 'H1(t)', 'Q1(t)', 'R1(t)', 'D1(t)')
% modello Gatto 3 popolazione (seconda dose)
x20 = 1e-6*ones(6,1);
time1 = 0:29;
tic
[t,x_vaccini3]= ode45('gatto_vaccini3', time, x20); 
toc
%
figure(5)
plot(t, x_vaccini3(:,:))
legend('S2(t)','E2(t)', 'P2(t)', 'I2(t)', 'A2(t)','R2(t)')

figure(6)
plot(t, x_vaccini3(:,2:end))
legend('E2(t)', 'P2(t)', 'I2(t)', 'A2(t)', 'R2(t)')

%% inclusione realismo con cascate di ODE su E ed H come a fine pagina 11 dell'articolo
global xc0 xc10 xc20
xc0 = [x0, zeros(1,3)];
tic
[t,x_vaccini]= ode45('gatto_vaccini_cascate', time, xc0); 
toc
%vettore di plot, non prendo gli stadi intermedi ma solo quell finali
plot_vett = [1 3 4 5 6 9 10 11 12];
figure(1)
plot(t, x_vaccini(:,plot_vett))
legend('S(t)','E(t)', 'P(t)', 'I(t)', 'A(t)', 'H(t)', 'Q(t)', 'R(t)', 'D(t)')
figure(2)
plot(t, x_vaccini(:,2:end))
legend('E(t)', 'P(t)', 'I(t)', 'A(t)', 'H(t)', 'Q(t)', 'R(t)', 'D(t)')
% modello Gatto 2 popolazione (prima dose)
xc10 = 1e-9*ones(12,1);
tic
[t,x_vaccini2]= ode45('gatto_vaccini_cascate2', time, xc10); 
toc
%
figure(3)
plot(t, x_vaccini2(:,:))
legend('S1(t)','E1(t)', 'P1(t)', 'I1(t)', 'A1(t)', 'H1(t)', 'Q1(t)', 'R1(t)', 'D1(t)')
figure(4)
plot(t, x_vaccini2(:,2:end))
legend('E1(t)', 'P1(t)', 'I1(t)', 'A1(t)', 'H1(t)', 'Q1(t)', 'R1(t)', 'D1(t)')
% modello Gatto 3 popolazione (seconda dose)
xc20 = 1e-6*ones(7,1);
time1 = 0:29;
tic
[t,x_vaccini3]= ode45('gatto_vaccini_cascate3', time, xc20); 
toc
%
figure(5)
plot(t, x_vaccini3(:,:))
legend('S2(t)','E2(t)', 'P2(t)', 'I2(t)', 'A2(t)','R2(t)')

figure(6)
plot(t, x_vaccini3(:,2:end))
legend('E2(t)', 'P2(t)', 'I2(t)', 'A2(t)', 'R2(t)')
%% VEDIAMO CHE SUCCEDE ACCOPPIANDO LE DINAMICHE
x0tot = [0.99 0.01 zeros(1,22)];
tic
[t,x_vaccini_tot]= ode45('gatto_vaccini_TOT', time, x0tot); 
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
legend('E2(t)', 'P2(t)', 'I2(t)', 'A2(t)','R2(t)')
%%