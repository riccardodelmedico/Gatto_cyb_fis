clc
clear
clear all
%%
parameters_vaccini;
dati_vaccini;
global lambda deltaE deltaP sigma eta gammaI alfaI gammaA zeta gammaH ...
       alfaH gammaQ gammaA Lvect x0 N eff1 eff2 ef1
 
format longg
N=136;
time= 0:1:N-1;
E0= 10/pop;
x0= [1-E0 E0 0 0 0 0 0 0 0]; 
%% modello Gatto con lockdown
tic
[t,x_vaccini]= ode45('gatto_vaccini', time, x0); 
toc

figure(2)
plot(t, x_vaccini(:,:))
legend('S(t)', 'E(t)', 'P(t)', 'I(t)', 'A(t)', 'H(t)', 'Q(t)', 'R(t)', 'D(t)')
