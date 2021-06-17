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

format longg
parameters_vaccini;
dati_vaccini;

Lvect = zeros(1,N);
x0= [1-1*E0 1*E0 zeros(1,22)];

tic
[x_vaccini_tot]= ode4(@gatto_vaccini_unico, 0, 1, 30, x0'); 
toc
% lambda in questo file accoppia le sottodinamiche dei 3 set di equazioni dei vaccini
soluzione= zeros(30, 24);
for j=1:1:30
    for i= 1:1:24
        soluzione(j,i)=x_vaccini_tot(24*(j-1)+i);
    end
end


E = [soluzione(:,2)];
P = [soluzione(:,3)];
I = [soluzione(:,4)];
A = [soluzione(:,5)];

figure(1)
plot(log(E))
hold on
plot(log(P))
plot(log(I))
plot(log(A))
hold off
legend('E', 'P','I', 'A')
xlabel('Days')
ylabel('Logarithmic scale')

figure(2)
plot(E)
hold on
plot(P)
plot(I)
plot(A)
hold off
legend('E', 'P','I', 'A')
xlabel('Days')

tasso_expo = log(P(end)) - log(P(end-1))
t_radd = log(2)/tasso_expo


%% ritaratura R0, mettendo il tempo di raddoppio a CIRCA 3, A ''OCCHIO''

parameters_vaccini_R0_raddoppio;


dati_vaccini;

Lvect = zeros(1,N);
x0= [1-1*E0 1*E0 zeros(1,22)];

tic
[x_vaccini_tot]= ode4(@gatto_vaccini_unico, 0, 1, 30, x0'); 
toc
% lambda in questo file accoppia le sottodinamiche dei 3 set di equazioni dei vaccini
soluzione= zeros(30, 24);
for j=1:1:30
    for i= 1:1:24
        soluzione(j,i)=x_vaccini_tot(24*(j-1)+i);
    end
end


E = [soluzione(:,2)];
P = [soluzione(:,3)];
I = [soluzione(:,4)];
A = [soluzione(:,5)];

figure(1)
plot(log(E))
hold on
plot(log(P))
plot(log(I))
plot(log(A))
hold off
legend('E', 'P','I', 'A')
xlabel('Days')
ylabel('Logarithmic scale')

figure(2)
plot(E)
hold on
plot(P)
plot(I)
plot(A)
hold off
legend('E', 'P','I', 'A')
xlabel('Days')

tasso_expo = log(P(end)) - log(P(end-1))
t_radd = log(2)/tasso_expo




%% rifacciamo l'analisi del tasso di raddoppio nel nuovo modello per vedere se c'è stato un rallentamento 
%causa cambio ipotes i sulle distribuzioni degli infetti (da esponenziale a gamma) 

%rifacciamo lo studio originale

parameters_vaccini

dati_vaccini;

Lvect = zeros(1,N);
x0= [1-1*E0 1*E0 zeros(1,22)];

tic
[x_vaccini_tot]= ode4(@gatto_vaccini_unico_cascate, 0, 1, 30, x0_casc'); 
toc
% lambda in questo file accoppia le sottodinamiche dei 3 set di equazioni dei vaccini
soluzione= zeros(30,31 );
for j=1:1:30
    for i= 1:1:31
        soluzione(j,i)=x_vaccini_tot(31*(j-1)+i);
    end
end


E = [soluzione(:,3)];
P = [soluzione(:,4)];
I = [soluzione(:,5)];
A = [soluzione(:,6)];

figure(1)
plot(log(E))
hold on
plot(log(P))
plot(log(I))
plot(log(A))
hold off
legend('E', 'P','I', 'A')
xlabel('Days')
ylabel('Logarithmic scale')

figure(2)
plot(E)
hold on
plot(P)
plot(I)
plot(A)
hold off
legend('E', 'P','I', 'A')
xlabel('Days')

tasso_expo = log(P(end)) - log(P(end-1))
t_radd = log(2)/tasso_expo

%% cerchiamo di tarare allora 
parameters_vaccini_R0_raddoppio
dati_vaccini;

Lvect = zeros(1,N);
x0= [1-1*E0 1*E0 zeros(1,22)];

tic
[x_vaccini_tot]= ode4(@gatto_vaccini_unico_cascate, 0, 1, 30, x0_casc'); 
toc
% lambda in questo file accoppia le sottodinamiche dei 3 set di equazioni dei vaccini
soluzione= zeros(30,31 );
for j=1:1:30
    for i= 1:1:31
        soluzione(j,i)=x_vaccini_tot(31*(j-1)+i);
    end
end


E = [soluzione(:,3)];
P = [soluzione(:,4)];
I = [soluzione(:,5)];
A = [soluzione(:,6)];

figure(1)
plot(log(E))
hold on
plot(log(P))
plot(log(I))
plot(log(A))
hold off
legend('E', 'P','I', 'A')
xlabel('Days')
ylabel('Logarithmic scale')

figure(2)
plot(E)
hold on
plot(P)
plot(I)
plot(A)
hold off
legend('E', 'P','I', 'A')
xlabel('Days')

tasso_expo = log(P(end)) - log(P(end-1))
t_radd = log(2)/tasso_expo



%% rifacciamo l'analisi del tasso di raddoppio nel nuovo modello per vedere se c'è stato un rallentamento 
%causa cambio ipotes i sulle distribuzioni degli infetti (da esponenziale a gamma) 
parameters_vaccini

dati_vaccini;

Lvect = zeros(1,N);
% x0_casc_inf= [1-1*E0 1*E0 zeros(1,22)];

tic
[x_vaccini_tot]= ode4(@gatto_vaccini_unico_cascatesoloInfetti, 0, 1, 30, x0_casc_inf'); 
toc
% lambda in questo file accoppia le sottodinamiche dei 3 set di equazioni dei vaccini
soluzione= zeros(30,42 );
for j=1:1:30
    for i= 1:1:42
        soluzione(j,i)=x_vaccini_tot(42*(j-1)+i);
    end
end


E = [soluzione(:,2)];
P = [soluzione(:,3)];
I = [soluzione(:,4)];
A = [soluzione(:,7)];

figure(1)
plot(log(E))
hold on
plot(log(P))
plot(log(I))
plot(log(A))
hold off
legend('E', 'P','I', 'A')
xlabel('Days')
ylabel('Logarithmic scale')

figure(2)
plot(E)
hold on
plot(P)
plot(I)
plot(A)
hold off
legend('E', 'P','I', 'A')
xlabel('Days')

tasso_expo = log(P(end)) - log(P(end-1))
t_radd = log(2)/tasso_expo

%% ritariamo anche questo

parameters_vaccini_R0_raddoppio


dati_vaccini;

Lvect = zeros(1,N);
% x0_casc_inf= [1-1*E0 1*E0 zeros(1,22)];

tic
[x_vaccini_tot]= ode4(@gatto_vaccini_unico_cascatesoloInfetti, 0, 1, 30, x0_casc_inf'); 
toc
% lambda in questo file accoppia le sottodinamiche dei 3 set di equazioni dei vaccini
soluzione= zeros(30,42 );
for j=1:1:30
    for i= 1:1:42
        soluzione(j,i)=x_vaccini_tot(42*(j-1)+i);
    end
end


E = [soluzione(:,2)];
P = [soluzione(:,3)];
I = [soluzione(:,4)];
A = [soluzione(:,7)];

figure(1)
plot(log(E))
hold on
plot(log(P))
plot(log(I))
plot(log(A))
hold off
legend('E', 'P','I', 'A')
xlabel('Days')
ylabel('Logarithmic scale')

figure(2)
plot(E)
hold on
plot(P)
plot(I)
plot(A)
hold off
legend('E', 'P','I', 'A')
xlabel('Days')

tasso_expo = log(P(end)) - log(P(end-1))
t_radd = log(2)/tasso_expo